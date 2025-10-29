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

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <arkode/arkode.hpp>
#include <arkode/arkode_splittingstep.hpp>
#include <sundials/sundials_core.hpp>

#include "arkode_mristep_impl.h"

#include "sundials4py_types.hpp"

namespace nb = nanobind;
using sundials::experimental::our_make_shared;
using sundials::experimental::SplittingStepCoefficientsDeleter;

namespace sundials4py {

void bind_arkode_splittingstep(nb::module_& m)
{
#include "arkode_splittingstep_generated.hpp"

  m.def("SplittingStepCreate",
        [](std::vector<SUNStepper> steppers, int partitions, sunrealtype t0,
           N_Vector y0, SUNContext sunctx) -> void* {
          return SplittingStepCreate(steppers.data(), partitions, t0, y0, sunctx);
        });

  m.def("SplittingStepReInit",
        [](void* arkode_mem, std::vector<SUNStepper> steppers, int partitions,
           sunrealtype t0, N_Vector y0) -> int
        {
          return SplittingStepReInit(arkode_mem, steppers.data(), partitions,
                                     t0, y0);
        });
}

} // namespace sundials4py
