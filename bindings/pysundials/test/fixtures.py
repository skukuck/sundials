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

from pysundials.core import *


@pytest.fixture
def sunctx():
    sunctx = SUNContextView.Create()
    yield sunctx


@pytest.fixture
def nvec(sunctx):
    nvec = NVectorView.Create(N_VNew_Serial(10, sunctx.get()))
    yield nvec


@pytest.fixture
def sunstepper(sunctx):
    status, stepper1 = SUNStepper_Create(sunctx.get())
    stepper = SUNStepperView.Create(stepper1)

    # Dummy callback implementations
    def evolve_fn(stepper, tout, vret, tret):
        return 0

    def one_step_fn(stepper, tout, vret, tret):
        return 0

    def full_rhs_fn(stepper, t, v, f, mode):
        return 0

    def reinit_fn(stepper, t, y):
        return 0

    def reset_fn(stepper, t, y):
        return 0

    def reset_ckpt_idx_fn(stepper, idx):
        return 0

    def stop_time_fn(stepper, tstop):
        return 0

    def step_direction_fn(stepper, direction):
        return 0

    def forcing_fn(stepper, tshift, tscale, forcing, nforcing):
        return 0

    def get_num_steps_fn(stepper, nst):
        return 0

    # Set all function pointers
    SUNStepper_SetEvolveFn(stepper.get(), evolve_fn)
    SUNStepper_SetOneStepFn(stepper.get(), one_step_fn)
    SUNStepper_SetFullRhsFn(stepper.get(), full_rhs_fn)
    SUNStepper_SetReInitFn(stepper.get(), reinit_fn)
    SUNStepper_SetResetFn(stepper.get(), reset_fn)
    SUNStepper_SetResetCheckpointIndexFn(stepper.get(), reset_ckpt_idx_fn)
    SUNStepper_SetStopTimeFn(stepper.get(), stop_time_fn)
    SUNStepper_SetStepDirectionFn(stepper.get(), step_direction_fn)
    SUNStepper_SetForcingFn(stepper.get(), forcing_fn)
    SUNStepper_SetGetNumStepsFn(stepper.get(), get_num_steps_fn)

    yield stepper
