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
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/tuple.h>

#include <sundials/sundials_core.hpp>
#include <sundials/sundials_adjointcheckpointscheme.hpp>
#include <sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h>

#include "sundials/sundials_adjointcheckpointscheme.h"
#include "sundials_adjointcheckpointscheme_impl.h"

#include "sundials4py_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunadjointcheckpointscheme_fixed(nb::module_& m)
{
  m.def(
    "SUNAdjointCheckpointScheme_Create_Fixed",
    [](SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, suncountertype interval,
       suncountertype estimate, sunbooleantype keep, SUNContext sunctx)
    {
      SUNAdjointCheckpointScheme check_scheme = nullptr;

      int status = SUNAdjointCheckpointScheme_Create_Fixed(io_mode, mem_helper,
                                                           interval, estimate,
                                                           keep, sunctx,
                                                           &check_scheme);

      return std::make_tuple(status, our_make_shared<std::remove_pointer_t<SUNAdjointCheckpointScheme>, SUNAdjointCheckpointSchemeDeleter>(check_scheme));
    },
    sundials4py::keep_alive_tuple<1, 5>());
}

} // namespace sundials4py
