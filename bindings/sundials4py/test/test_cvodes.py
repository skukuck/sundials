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
import numpy as np
from fixtures import *
from sundials4py.core import *
from sundials4py.cvodes import *
from problems import AnalyticODE


@pytest.mark.skipif(
    sunrealtype == np.float32, reason="Test not supported for sunrealtype=np.float32"
)
def test_cvodes_ivp(sunctx):
    y = N_VNew_Serial(1, sunctx)
    ls = SUNLinSol_SPGMR(y, 0, 0, sunctx)

    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y)

    cvode = CVodeCreate(CV_BDF, sunctx)

    status = CVodeInit(cvode.get(), lambda t, y, ydot, _: ode_problem.f(t, y, ydot), 0, y)
    assert status == CV_SUCCESS

    status = CVodeSStolerances(cvode.get(), SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == CV_SUCCESS

    status = CVodeSetLinearSolver(cvode.get(), ls, None)
    assert status == CV_SUCCESS

    tout = 10.0
    status, tret = CVode(cvode.get(), tout, y, CV_NORMAL)
    assert status == CV_SUCCESS

    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)
    assert np.allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=1e-2)

    status, num_steps = CVodeGetNumSteps(cvode.get())
    assert status == CV_SUCCESS
    assert num_steps > 0


def test_cvodes_fsa(sunctx):
    # TODO(CJB): implement
    pass


def test_cvodes_adjoint(sunctx):
    # TODO(CJB): implement
    pass
