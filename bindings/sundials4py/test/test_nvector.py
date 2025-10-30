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


@pytest.mark.parametrize("vector_type", ["serial"])
def test_create_nvector(vector_type, sunctx):
    if vector_type == "serial":
        nvec = N_VNew_Serial(5, sunctx)
    else:
        raise ValueError("Unknown vector type")
    assert nvec is not None

    arr = N_VGetArrayPointer(nvec)
    assert arr.shape[0] == 5

    arr[:] = np.array([5.0, 4.0, 3.0, 2.0, 1.0], dtype=sunrealtype)
    assert np.allclose(N_VGetArrayPointer(nvec), [5.0, 4.0, 3.0, 2.0, 1.0])

    N_VConst(2.0, nvec)
    assert np.allclose(arr, 2.0)


@pytest.mark.parametrize("vector_type", ["serial"])
def test_make_nvector(vector_type, sunctx):
    arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=sunrealtype)
    if vector_type == "serial":
        nvec = N_VMake_Serial(5, arr, sunctx)
    else:
        raise ValueError("Unknown vector type")
    assert nvec is not None

    N_VConst(2.0, nvec)
    assert np.allclose(arr, 2.0)

    arr[:] = np.array([5.0, 4.0, 3.0, 2.0, 1.0], dtype=sunrealtype)
    assert np.allclose(N_VGetArrayPointer(nvec), [5.0, 4.0, 3.0, 2.0, 1.0])


@pytest.mark.parametrize("vector_type", ["serial"])
def test_setarraypointer(vector_type, sunctx):
    if vector_type == "serial":
        nvec = N_VNew_Serial(5, sunctx)
    else:
        raise ValueError("Unknown vector type")
    assert nvec is not None

    arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=sunrealtype)
    N_VSetArrayPointer(arr, nvec)

    assert np.allclose(N_VGetArrayPointer(nvec), arr)

    N_VScale(2.0, nvec, nvec)

    assert np.allclose(arr, [2.0, 4.0, 6.0, 8.0, 10.0])


# Test an operation that involves vector arrays
@pytest.mark.parametrize("vector_type", ["serial"])
def test_nvlinearcombination(vector_type, sunctx):
    if vector_type == "serial":
        nvec1 = N_VNew_Serial(5, sunctx)
        nvec2 = N_VNew_Serial(5, sunctx)
    else:
        raise ValueError("Unknown vector type")

    arr1 = N_VGetArrayPointer(nvec1)
    arr1[:] = [1.0, 2.0, 3.0, 4.0, 5.0]

    arr2 = N_VGetArrayPointer(nvec2)
    arr2[:] = [10.0, 20.0, 30.0, 40.0, 50.0]

    c = np.array([1.0, 0.1], dtype=sunrealtype)
    X = [nvec1, nvec2]

    z = N_VNew_Serial(5, sunctx)
    N_VConst(0.0, z)

    N_VLinearCombination(2, c, X, z)

    assert np.allclose(N_VGetArrayPointer(z), [2.0, 4.0, 6.0, 8.0, 10.0])
