/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *------------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2025, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_TYPES_HPP
#define _SUNDIALS4PY_TYPES_HPP

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

namespace sundials4py {

using Array1d = nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig>;
using CpuArray1d =
  nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

} // namespace sundials4py

#endif
