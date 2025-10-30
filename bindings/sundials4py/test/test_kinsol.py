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
from sundials4py.kinsol import *
from problems import AnalyticNonlinearSys


def test_kinsol(sunctx):
    NEQ = 3
    m_aa = 2
    tol = 1e-4
    kin_view = KINCreate(sunctx)
    problem = AnalyticNonlinearSys(None)
    u = N_VNew_Serial(NEQ, sunctx)

    def fp_function(u, g, _):
        return problem.fixed_point_fn(u, g)

    kin_status = KINSetMAA(kin_view.get(), m_aa)
    assert kin_status == KIN_SUCCESS
    kin_status = KINInit(kin_view.get(), fp_function, u)
    assert kin_status == KIN_SUCCESS
    kin_status = KINSetFuncNormTol(kin_view.get(), tol)
    assert kin_status == KIN_SUCCESS

    # initial guess
    u_data = N_VGetArrayPointer(u)
    u_data[:] = [0.1, 0.1, -0.1]

    # no scaling used
    scale = N_VNew_Serial(NEQ, sunctx)
    N_VConst(1.0, scale)

    kin_status = KINSol(kin_view.get(), u, KIN_FP, scale, scale)
    assert kin_status == KIN_SUCCESS

    u_expected = N_VNew_Serial(NEQ, sunctx)
    u_expected_data = N_VGetArrayPointer(u_expected)
    problem.solution(u_expected)
    assert np.allclose(u_data, u_expected_data, atol=1e-6)
