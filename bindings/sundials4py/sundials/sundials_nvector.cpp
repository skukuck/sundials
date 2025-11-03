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
 * SUNDIALS N_Vector class. It contains hand-written code for functions
 * that require special treatment, and includes the generated code
 * produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_nvector.hpp>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_nvector(nb::module_& m)
{
#include "sundials_nvector_generated.hpp"

  m.def("N_VGetArrayPointer",
        [](N_Vector v)
        {
          auto ptr = N_VGetArrayPointer(v);
          if (!ptr)
          {
            throw sundials4py::error_returned("Failed to get array pointer");
          }
          auto owner = nb::find(v);
          size_t shape[1]{static_cast<size_t>(N_VGetLength(v))};
          return sundials4py::Array1d(ptr, 1, shape, owner);
        }, nb::rv_policy::reference);

  m.def("N_VGetDeviceArrayPointer",
        [](N_Vector v)
        {
          auto ptr = N_VGetDeviceArrayPointer(v);
          if (!ptr)
          {
            throw sundials4py::error_returned("Failed to get array pointer");
          }
          auto owner = nb::find(v);
          size_t shape[1]{static_cast<size_t>(N_VGetLength(v))};
          return sundials4py::Array1d(ptr, 1, shape, owner);
        }, nb::rv_policy::reference);

  m.def("N_VSetArrayPointer",
        [](sundials4py::Array1d arr, N_Vector v)
        {
          if (arr.shape(0) != N_VGetLength(v))
          {
            throw sundials4py::error_returned(
              "Array shape does not match vector length");
          }
          N_VSetArrayPointer(arr.data(), v);
        });
}

} // namespace sundials4py
