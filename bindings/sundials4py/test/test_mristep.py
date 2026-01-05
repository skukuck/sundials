#!/bin/python
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2025, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------


import pytest
import weakref
import numpy as np
from numpy.testing import assert_allclose
from sundials4py.core import *
from sundials4py.arkode import *
from problems import AnalyticMultiscaleODE
from fixtures import sunctx


def test_multirate(sunctx):
    ode_problem = AnalyticMultiscaleODE()

    t0, tf = AnalyticMultiscaleODE.T0, 0.01

    def fslow(t, y, ydot, _):
        return ode_problem.f_linear(t, y, ydot, None)

    def ffast(t, y, ydot, _):
        # # TODO(CJB): fix MRIStepInnerStepper_GetForcingData
        # inner_stepper = ode_problem.inner_stepper
        # # test MRIStepInnerStepper_GetForcingData
        # status, tshift, tscale, forcing, nforcing = MRIStepInnerStepper_GetForcingData(inner_stepper)
        # assert status == ARK_SUCCESS
        # assert len(forcing) == nforcing

        return ode_problem.f_nonlinear(t, y, ydot, None)

    y = N_VNew_Serial(1, sunctx)
    y0 = N_VClone(y)

    ode_problem.set_init_cond(y)
    ode_problem.set_init_cond(y0)

    # create fast integrator
    inner_ark = ERKStepCreate(ffast, t0, y, sunctx)
    status = ARKodeSetFixedStep(inner_ark.get(), 5e-3)
    assert status == ARK_SUCCESS

    status, inner_stepper = ARKodeCreateMRIStepInnerStepper(inner_ark.get())
    assert status == ARK_SUCCESS

    # store inner_stepper in ode_problem so we can access it in ffast
    ode_problem.inner_stepper = inner_stepper

    # create slow integrator
    ark = MRIStepCreate(fslow, None, t0, y, inner_stepper, sunctx)
    status = ARKodeSetFixedStep(ark.get(), 1e-3)
    assert status == ARK_SUCCESS

    tout = tf
    status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
    assert status == ARK_SUCCESS

    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)
    assert_allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=1e-2)

    # We must set this to None to ensure inner_stepper can be garbage collected
    # If we do not do this, then nanobind will warn that references are leaked.
    # This seems to be unavoidable without setting this to None or using a weakref.
    # Its possible newer versions of Python may not result in the warning.
    ode_problem.inner_stepper = None
