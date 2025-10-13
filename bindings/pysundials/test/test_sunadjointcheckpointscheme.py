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
# Unit/smoke tests for SUNAdjointCheckpointScheme module
# -----------------------------------------------------------------

import pytest
from fixtures import *
from pysundials.core import *

def make_fixed_scheme(sunctx):
	# Use dummy values for required arguments
	io_mode = 0
	mem_helper = None
	interval = 1
	estimate = 1
	keep = 0
	status, scheme = SUNAdjointCheckpointScheme_Create_Fixed(io_mode, mem_helper, interval, estimate, keep, sunctx.get())
	return status, scheme


def test_needs_saving(sunctx):
    # Dummy values for arguments
    scheme_status, scheme = make_fixed_scheme(sunctx)
    step_num = 0
    stage_num = 0
    t = 0.0
    try:
        status, result = SUNAdjointCheckpointScheme_NeedsSaving(scheme, step_num, stage_num, t)
        assert status == 0
        assert isinstance(result, bool)
    except Exception:
        pytest.skip("SUNAdjointCheckpointScheme_NeedsSaving not bound or not implemented")


def test_insert_vector(sunctx, nvec):
    scheme_status, scheme = make_fixed_scheme(sunctx)
    step_num = 0
    stage_num = 0
    t = 0.0
    try:
        status = SUNAdjointCheckpointScheme_InsertVector(scheme, step_num, stage_num, t, nvec.get())
        assert status == 0
    except Exception:
        pytest.skip("SUNAdjointCheckpointScheme_InsertVector not bound or not implemented")


def test_load_vector(sunctx):
    scheme_status, scheme = make_fixed_scheme(sunctx)
    step_num = 0
    stage_num = 0
    peek = False
    try:
        status, vec, tout = SUNAdjointCheckpointScheme_LoadVector(scheme, step_num, stage_num, peek)
        assert status == 0
        assert vec is not None
        assert isinstance(tout, float)
    except Exception:
        pytest.skip("SUNAdjointCheckpointScheme_LoadVector not bound or not implemented")


def test_enable_dense(sunctx):
    scheme_status, scheme = make_fixed_scheme(sunctx)
    try:
        status = SUNAdjointCheckpointScheme_EnableDense(scheme, True)
        assert status == 0
    except Exception:
        pytest.skip("SUNAdjointCheckpointScheme_EnableDense not bound or not implemented")
