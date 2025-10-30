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
from fixtures import *
from sundials4py.core import *

# Note: some of these tests will fail if SUNDIALS error checks are turned on because
# we dont properly mock some of the requirements


def test_smoke_create_dense(sunctx, nvec):
    A = SUNDenseMatrix(2, 2, sunctx)
    LS = SUNLinSol_Dense(nvec, A, sunctx)
    assert LS is not None


def test_smoke_create_band(sunctx, nvec):
    A = SUNBandMatrix(2, 1, 1, sunctx)
    LS = SUNLinSol_Band(nvec, A, sunctx)
    assert LS is not None


def test_smoke_create_spgmr(sunctx, nvec):
    LS = SUNLinSol_SPGMR(nvec, SUN_PREC_NONE, 0, sunctx)
    assert LS is not None


def test_smoke_create_pcg(sunctx, nvec):
    LS = SUNLinSol_PCG(nvec, SUN_PREC_NONE, 0, sunctx)
    assert LS is not None


def test_smoke_create_spbcgs(sunctx, nvec):
    LS = SUNLinSol_SPBCGS(nvec, SUN_PREC_NONE, 0, sunctx)
    assert LS is not None


def test_smoke_create_sptfqmr(sunctx, nvec):
    LS = SUNLinSol_SPTFQMR(nvec, SUN_PREC_NONE, 0, sunctx)
    assert LS is not None


def test_smoke_get_type_and_id(sunctx, nvec):
    A = SUNDenseMatrix(2, 2, sunctx)
    LS = SUNLinSol_Dense(nvec, A, sunctx)
    typ = SUNLinSolGetType(LS)
    id_ = SUNLinSolGetID(LS)
    assert isinstance(typ, int)
    assert isinstance(id_, int)


def test_smoke_initialize_setup(sunctx, nvec):
    A = SUNDenseMatrix(2, 2, sunctx)
    LS = SUNLinSol_Dense(nvec, A, sunctx)
    ret_init = SUNLinSolInitialize(LS)
    ret_setup = SUNLinSolSetup(LS, A)
    assert isinstance(ret_init, int)
    assert isinstance(ret_setup, int)


def test_smoke_num_iters_resnorm_lastflag(sunctx, nvec):
    LS = SUNLinSol_SPGMR(nvec, 0, 0, sunctx)
    niters = SUNLinSolNumIters(LS)
    resnorm = SUNLinSolResNorm(LS)
    lastflag = SUNLinSolLastFlag(LS)
    assert isinstance(niters, int)
    assert isinstance(resnorm, float)
    assert isinstance(lastflag, int)
