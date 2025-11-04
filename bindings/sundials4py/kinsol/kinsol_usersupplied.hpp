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

#ifndef _SUNDIALS4PY_KINSOL_USERSUPPLIED_HPP
#define _SUNDIALS4PY_KINSOL_USERSUPPLIED_HPP

#include <cstdlib>
#include <cstring>

#include "sundials/sundials_types.h"
#include "sundials4py.hpp"

#include <kinsol/kinsol.h>
#include <kinsol/kinsol_ls.h>

#include <sundials/sundials_core.hpp>

#include "sundials4py_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

///////////////////////////////////////////////////////////////////////////////
// KINSOL user-supplied function table
// Every package-level user-supplied function must be in this table.
// The user-supplied function table is passed to KINSOL as user_data.
///////////////////////////////////////////////////////////////////////////////

struct kinsol_user_supplied_fn_table
{
  // KINSOL user-supplied function pointers
  nb::object sysfn, dampingfn, depthfn;

  // KINSOL LS user-supplied function pointers
  nb::object lsjacfn, lsjactimesvecfn, lsjtvsysfn, lsprecsetupfn, lsprecsolvefn;
};

///////////////////////////////////////////////////////////////////////////////
// KINSOL user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline kinsol_user_supplied_fn_table* kinsol_user_supplied_fn_table_alloc()
{
  // We must use malloc since KINFree calls free
  auto fn_table = static_cast<kinsol_user_supplied_fn_table*>(
    std::malloc(sizeof(kinsol_user_supplied_fn_table)));

  // Zero out the memory
  std::memset(fn_table, 0, sizeof(kinsol_user_supplied_fn_table));

  return fn_table;
}

template<typename... Args>
inline int kinsol_sysfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<KINSysFn>, kinsol_user_supplied_fn_table,
    1>(&kinsol_user_supplied_fn_table::sysfn, std::forward<Args>(args)...);
}

using KINDampingStdFn = std::tuple<int, sunrealtype>(
  long int iter, N_Vector u_val, N_Vector g_val, sundials4py::Array1d qt_fn,
  long int depth, void* user_data);

inline int kinsol_dampingfn_wrapper(long int iter, N_Vector u_val, N_Vector g_val,
                                    sunrealtype* qt_fn_1d, long int depth,
                                    void* user_data, sunrealtype* damping_factor)
{
  auto fn_table = static_cast<kinsol_user_supplied_fn_table*>(user_data);
  auto fn       = nb::cast<std::function<KINDampingStdFn>>(fn_table->dampingfn);

  sundials4py::Array1d qt_fn(qt_fn_1d, {static_cast<unsigned long>(depth)});

  auto result = fn(iter, u_val, g_val, qt_fn, depth, nullptr);

  *damping_factor = std::get<1>(result);

  return std::get<0>(result);
}

using KINDepthStdFn = std::tuple<int, long int>(
  long int iter, N_Vector u_val, N_Vector g_val, N_Vector f_val,
  std::vector<N_Vector> df, sundials4py::Array1d R_mat, long int depth,
  void* user_data, std::vector<sunbooleantype> remove_indices);

inline int kinsol_depthfn_wrapper(long int iter, N_Vector u_val, N_Vector g_val,
                                  N_Vector f_val, N_Vector* df_1d,
                                  sunrealtype* R_mat_1d, long int depth,
                                  void* user_data, long int* new_depth,
                                  sunbooleantype* remove_indices_1d)
{
  auto fn_table = static_cast<kinsol_user_supplied_fn_table*>(user_data);
  auto fn       = nb::cast<std::function<KINDepthStdFn>>(fn_table->depthfn);

  std::vector<N_Vector> df(df_1d, df_1d + depth);
  sundials4py::Array1d R_mat(R_mat_1d,
                             {static_cast<unsigned long>(depth * depth)});
  // TODO(CJB): enable this if it becomes used in the future
  std::vector<sunbooleantype> remove_indices(0);

  auto result = fn(iter, u_val, g_val, f_val, df, R_mat, depth, nullptr,
                   remove_indices);

  *new_depth = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int kinsol_lsjacfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<KINLsJacFn>, kinsol_user_supplied_fn_table,
    3>(&kinsol_user_supplied_fn_table::lsjacfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int kinsol_lsprecsetupfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<KINLsPrecSetupFn>, kinsol_user_supplied_fn_table,
    1>(&kinsol_user_supplied_fn_table::lsprecsetupfn,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int kinsol_lsprecsolvefn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<KINLsPrecSolveFn>, kinsol_user_supplied_fn_table,
    1>(&kinsol_user_supplied_fn_table::lsprecsolvefn,
       std::forward<Args>(args)...);
}

using KINLsJacTimesVecStdFn = std::tuple<int, sunbooleantype>(N_Vector v,
                                                              N_Vector Jv,
                                                              N_Vector u,
                                                              void* user_data);

inline int kinsol_lsjactimesvecfn_wrapper(N_Vector v, N_Vector Jv, N_Vector u,
                                          sunbooleantype* new_u, void* user_data)
{
  auto fn_table = static_cast<kinsol_user_supplied_fn_table*>(user_data);
  auto fn =
    nb::cast<std::function<KINLsJacTimesVecStdFn>>(fn_table->lsjactimesvecfn);

  auto result = fn(v, Jv, u, nullptr);

  *new_u = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int kinsol_lsjtvsysfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<KINSysFn>, kinsol_user_supplied_fn_table,
    1>(&kinsol_user_supplied_fn_table::lsjtvsysfn, std::forward<Args>(args)...);
}

#endif // _SUNDIALS4PY_KINSOL_USERSUPPLIED_HPP
