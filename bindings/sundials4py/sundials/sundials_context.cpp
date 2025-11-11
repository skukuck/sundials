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
 * SUNDIALS SUNContext class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_context.hpp>
#include <sundials/sundials_logger.hpp>
#include <sundials/sundials_profiler.hpp>

#include "sundials/sundials_types.h"
#include "sundials_logger_impl.h"
#include "sundials_profiler_impl.h"

namespace nb = nanobind;
using namespace sundials::experimental;
using sundials::SUNContextDeleter;

#include "sundials_context_usersupplied.hpp"

namespace sundials4py {

using namespace sundials::experimental;

void bind_suncontext(nb::module_& m)
{
#include "sundials_context_generated.hpp"

  nb::class_<SUNContext_>(m, "SUNContext_");

  // Note: only one error handler can be pushed from python
  m.def("SUNContext_PushErrHandler",
        [](SUNContext sunctx,
           std::function<std::remove_pointer_t<SUNErrHandlerFn>> err_fn)
        {
          if (!sunctx->python)
          {
            sunctx->python = SUNContextFunctionTable_Alloc();
          }
          auto fn_table = static_cast<SUNContextFunctionTable*>(sunctx->python);

          if (fn_table->err_handler)
          {
            throw sundials4py::illegal_value(
              "SUNContext_PushErrHandler was already called. sundials4py only "
              "allows for SUNContext_PushErrHandler to be called once. Call "
              "SUNContext_PopErrHandler first, then call "
              "SUNContext_PushErrHandler again.")
          }

          fn_table->err_handler = nb::cast(err_fn);

          return SUNContext_PushErrHandler(sunctx, suncontext_errhandler_wrapper,
                                           sunctx->python);
        });
}

} // namespace sundials4py
