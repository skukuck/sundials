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

#ifndef _SUNDIALS4PY_SUNCONTEXT_USERSUPPLIED_HPP
#define _SUNDIALS4PY_SUNCONTEXT_USERSUPPLIED_HPP

#include <cstdlib>
#include <cstring>
#include <vector>

#include <sundials/sundials_context.h>
#include "sundials4py.hpp"
#include "sundials4py_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

// Function table for user-supplied error handler for SUNContext
struct SUNContextFunctionTable
{
  nb::object err_handler;
};

inline SUNContextFunctionTable* SUNContextFunctionTable_Alloc()
{
  auto fn_table = static_cast<SUNContextFunctionTable*>(
    std::malloc(sizeof(SUNContextFunctionTable)));
  std::memset(fn_table, 0, sizeof(SUNContextFunctionTable));
  return fn_table;
}

inline void suncontext_errhandler_wrapper(int line, const char* func,
                                          const char* file, const char* msg,
                                          SUNErrCode err_code,
                                          void* err_user_data, SUNContext sunctx)
{
  auto fn_table = static_cast<SUNContextFunctionTable*>(err_user_data);
  fn_table->err_handler(line, func, file, msg, err_code, nullptr, sunctx);
}

#endif