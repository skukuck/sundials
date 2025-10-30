/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This file is the entrypoint for the Python binding code for the
 * SUNDIALS SUNAdjointCheckpointScheme class. It contains hand-written
 * code for functions that require special treatment, and includes the
 * generated code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_adjointcheckpointscheme.hpp>
#include <sundials/sundials_nvector.hpp>

#include "sundials_adjointcheckpointscheme_impl.h"

namespace nb = nanobind;
using namespace sundials::experimental;

using namespace sundials::experimental;

namespace sundials4py {

void bind_sunadjointcheckpointscheme(nb::module_& m)
{
#include "sundials_adjointcheckpointscheme_generated.hpp"

  nb::class_<SUNAdjointCheckpointScheme_>(m, "SUNAdjointCheckpointScheme_");

  nb::class_<SUNAdjointCheckpointSchemeView>(m,
                                             "SUNAdjointCheckpointSchemeView")
    .def_static("Create",
                &SUNAdjointCheckpointSchemeView::Create<SUNAdjointCheckpointScheme>)
    .def("get",
         nb::overload_cast<>(&SUNAdjointCheckpointSchemeView::get, nb::const_),
         nb::rv_policy::reference);
}

} // namespace sundials4py
