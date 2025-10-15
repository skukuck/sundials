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
# Unit/smoke tests for SUNAdjointStepper module
# -----------------------------------------------------------------

import pytest
from fixtures import *
from pysundials.core import *

def make_adjoint_stepper(sunctx, nvec):
    mem_helper = SUNMemoryHelperView.Create(SUNMemoryHelper_Sys(sunctx.get()))
    status, scheme = SUNAdjointCheckpointScheme_Create_Fixed(0, mem_helper.get(), 1, 1, 0, sunctx.get())
    scheme = SUNAdjointCheckpointSchemeView.Create(scheme)
    status, stepper1 = SUNStepper_Create(sunctx.get())
    status, stepper2 = SUNStepper_Create(sunctx.get())
    stepper1 = SUNStepperView.Create(stepper1)
    stepper2 = SUNStepperView.Create(stepper2)
    b1 = 0
    b2 = 0
    nsteps = 1
    t0 = 0.0
    y0 = nvec.get()
    adj_view = SUNAdjointStepperView.Create(stepper1.get(), b1, stepper2.get(), b2, nsteps, t0, y0, scheme.get(), sunctx.get())
    return adj_view, scheme, mem_helper, stepper1, stepper2

def test_create_adjoint_stepper(sunctx, nvec):
    adj_view, scheme, mem_helper, s1, s2 = make_adjoint_stepper(sunctx, nvec)
    assert adj_view.get() is not None

def test_adjointstepper_reinit(sunctx, nvec):
    adj_view, scheme, mem_helper, s1, s2 = make_adjoint_stepper(sunctx, nvec)
    t0 = 0.0
    y0 = nvec.get()
    tf = 1.0
    err = SUNAdjointStepper_ReInit(adj_view.get(), t0, y0, tf, nvec.get())
    assert isinstance(err, int)

def test_adjointstepper_evolve(sunctx, nvec):
    adj_view, scheme, mem_helper, s1, s2 = make_adjoint_stepper(sunctx, nvec)
    tout = 1.0
    sens = nvec.get()
    err, tret = SUNAdjointStepper_Evolve(adj_view.get(), tout, sens)
    assert isinstance(err, int)
    assert isinstance(tret, float)

def test_adjointstepper_onestep(sunctx, nvec):
    adj_view, scheme, mem_helper, s1, s2 = make_adjoint_stepper(sunctx, nvec)
    tout = 1.0
    sens = nvec.get()
    err, tret = SUNAdjointStepper_OneStep(adj_view.get(), tout, sens)
    assert isinstance(err, int)
    assert isinstance(tret, float)

def test_adjointstepper_recomputefwd(sunctx, nvec):
	adj_view, scheme, mem_helper, s1, s2 = make_adjoint_stepper(sunctx, nvec)
	start_idx = 0
	t0 = 0.0
	y0 = nvec.get()
	tf = 1.0
	err = SUNAdjointStepper_RecomputeFwd(adj_view.get(), start_idx, t0, y0, tf)
	assert isinstance(err, int)

# def test_adjointstepper_getnumsteps(sunctx, nvec):
# 	adj_view, scheme, mem_helper, s1, s2 = make_adjoint_stepper(sunctx, nvec)
# 	err, num_steps = SUNAdjointStepper_GetNumSteps(adj_view.get())
# 	assert isinstance(err, int)
# 	assert isinstance(num_steps, int)

# def test_adjointstepper_getnumrecompute(sunctx, nvec):
# 	adj_view, scheme, mem_helper, s1, s2 = make_adjoint_stepper(sunctx, nvec)
# 	err, num_recompute = SUNAdjointStepper_GetNumRecompute(adj_view.get())
# 	assert isinstance(err, int)
# 	assert isinstance(num_recompute, int)
