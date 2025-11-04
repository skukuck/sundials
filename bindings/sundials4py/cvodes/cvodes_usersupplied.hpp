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

#ifndef _SUNDIALS4PY_CVODE_USERSUPPLIED_HPP
#define _SUNDIALS4PY_CVODE_USERSUPPLIED_HPP

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_ls.h>

#include <sundials/sundials_core.hpp>

#include "sundials4py_helpers.hpp"

///////////////////////////////////////////////////////////////////////////////
// CVODE user-supplied function table
// Every integrator-level user-supplied function must be in this table.
// The user-supplied function table is passed to CVODE as user_data.
///////////////////////////////////////////////////////////////////////////////

struct cvode_user_supplied_fn_table
{
  // user-supplied function pointers
  nb::object f, rootfn, ewtn, rwtn, fNLS, projfn;

  // cvode_ls user-supplied function pointers
  nb::object lsjacfn, lsprecsetupfn, lsprecsolvefn, lsjactimessetupfn,
    lsjactimesvecfn, lslinsysfn, lsjacrhsfn;

  // cvode quadrature user-supplied function pointers
  nanobind::object fQ, fQS;

  // cvode FSA user-supplied function pointers
  nanobind::object fS, fS1;
};

///////////////////////////////////////////////////////////////////////////////
// CVODE user-supplied functions
///////////////////////////////////////////////////////////////////////////////

inline cvode_user_supplied_fn_table* cvode_user_supplied_fn_table_alloc()
{
  // We must use malloc since CVODEFree calls free
  auto fn_table = static_cast<cvode_user_supplied_fn_table*>(
    std::malloc(sizeof(cvode_user_supplied_fn_table)));

  // Zero out the memory
  std::memset(fn_table, 0, sizeof(cvode_user_supplied_fn_table));

  return fn_table;
}

template<typename... Args>
inline int cvode_f_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<std::remove_pointer_t<CVRhsFn>,
                                              cvode_user_supplied_fn_table,
                                              1>(&cvode_user_supplied_fn_table::f,
                                                 std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_rootfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVRootFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::rootfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_ewtfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVEwtFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::ewtn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_nlsrhsfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::fNLS, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjacfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacFn>, cvode_user_supplied_fn_table,
    4>(&cvode_user_supplied_fn_table::lsjacfn, std::forward<Args>(args)...);
}

using CVLsPrecSetupStdFn = std::tuple<int, sunbooleantype>(sunrealtype t, N_Vector y, 
      N_Vector fy, sunbooleantype jok,
      sunrealtype gamma, void *user_data);

inline int cvode_lsprecsetupfn_wrapper(sunrealtype t, N_Vector y, 
      N_Vector fy, sunbooleantype jok, sunbooleantype *jcurPtr, 
      sunrealtype gamma, void *user_data)
{
  auto fn_table = static_cast<cvode_user_supplied_fn_table*>(user_data);
  auto fn       = nb::cast<std::function<CVLsPrecSetupStdFn>>(fn_table->lsprecsetupfn);

  auto result = fn(t, y, fy, jok, gamma, nullptr);

  *jcurPtr = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int cvode_lsprecsolvefn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSolveFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsprecsolvefn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjactimessetupfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesSetupFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsjactimessetupfn,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjactimesvecfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesVecFn>, cvode_user_supplied_fn_table,
    2>(&cvode_user_supplied_fn_table::lsjactimesvecfn,
       std::forward<Args>(args)...);
}

using CVLsLinSysStdFn = std::tuple<int, sunbooleantype>(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix M, 
                                    sunbooleantype jok, sunrealtype gamma, 
                                    void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


inline int cvode_lslinsysfn_wrapper(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix M, 
                                    sunbooleantype jok, sunbooleantype *jcur, sunrealtype gamma, 
                                    void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  auto fn_table = static_cast<cvode_user_supplied_fn_table*>(user_data);
  auto fn       = nb::cast<std::function<CVLsLinSysStdFn>>(fn_table->lslinsysfn);

  auto result = fn(t, y, fy, M, jok, gamma, nullptr, tmp1, tmp2, tmp3);

  *jcur = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int cvode_lsjacrhsfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::lsjacrhsfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_projfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVProjFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::projfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_fQ_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVQuadRhsFn>, cvode_user_supplied_fn_table,
    1>(&cvode_user_supplied_fn_table::fQ, std::forward<Args>(args)...);
}

using CVQuadSensRhsStdFn = int(sunrealtype t, N_Vector y,
                               std::vector<N_Vector> yQS, N_Vector ydot,
                               std::vector<N_Vector> yQSdot, void* user_data,
                               N_Vector tmp, N_Vector tmpQ);

template<typename... Args>
inline int cvode_fQS_wrapper(Args... args)
{
  // return sundials4py::user_supplied_fn_caller<
  //   std::remove_pointer_t<CVQuadSensRhsFn>, cvode_user_supplied_fn_table,
  //   3>(&cvode_user_supplied_fn_table::fQS, std::forward<Args>(args)...);
}

using CVSensRhsStdFn = int(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                           std::vector<N_Vector> yS, std::vector<N_Vector> ySdot,
                           void* user_data, N_Vector tmp1, N_Vector tmp2);

inline int cvode_fS_wrapper(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                            N_Vector* yS, N_Vector* ySdot, void* user_data,
                            N_Vector tmp1, N_Vector tmp2)
{
  auto fn_table = static_cast<cvode_user_supplied_fn_table*>(user_data);
  auto fn       = nb::cast<std::function<CVSensRhsStdFn>>(fn_table->fS);

  std::vector<N_Vector> yS_1d(yS, yS + Ns);
  std::vector<N_Vector> ySdot_1d(ySdot, ySdot + Ns);

  return fn(Ns, t, y, ydot, yS_1d, ySdot_1d, nullptr, tmp1, tmp2);
}

template<typename... Args>
inline int cvode_fS1_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVSensRhs1Fn>, cvode_user_supplied_fn_table,
    3>(&cvode_user_supplied_fn_table::fS1, std::forward<Args>(args)...);
}

///////////////////////////////////////////////////////////////////////////////
// CVODE Adjoint user-supplied functions
///////////////////////////////////////////////////////////////////////////////

struct cvodea_user_supplied_fn_table
{
  // cvode adjoint user-supplied function pointers
  nb::object fB, fQB, fQBS;

  // cvode_ls adjoint user-supplied function pointers
  nb::object lsjacfnB, lsjacfnBS, lsprecsetupfnB, lsprecsetupfnBS,
    lsprecsolvefnB, lsprecsolvefnBS, lsjactimessetupfnB, lsjactimessetupfnBS,
    lsjactimesvecfnB, lsjactimesvecfnBS, lslinsysfnB, lslinsysfnBS;
};

inline cvodea_user_supplied_fn_table* cvodea_user_supplied_fn_table_alloc()
{
  // We must use malloc since CVODEFree calls free
  auto fn_table = static_cast<cvodea_user_supplied_fn_table*>(
    std::malloc(sizeof(cvodea_user_supplied_fn_table)));

  // Zero out the memory
  std::memset(fn_table, 0, sizeof(cvodea_user_supplied_fn_table));

  return fn_table;
}

template<typename... Args>
inline int cvode_fB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVRhsFnB>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::fB, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_fQB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVQuadRhsFnB>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::fQB, std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjacfnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacFnB>, cvodea_user_supplied_fn_table,
    4>(&cvodea_user_supplied_fn_table::lsjacfnB, std::forward<Args>(args)...);
}

using CVLsPrecSetupStdFnB = std::tuple<int, sunbooleantype>(sunrealtype t, N_Vector y,
   N_Vector yB, N_Vector fyB, sunbooleantype jokB,
   sunrealtype gammaB, void *user_dataB);

inline int cvode_lsprecsetupfnB_wrapper(sunrealtype t, N_Vector y,
   N_Vector yB, N_Vector fyB, sunbooleantype jokB, sunbooleantype *jcurPtrB,
   sunrealtype gammaB, void *user_dataB)
{
  auto fn_table = static_cast<cvodea_user_supplied_fn_table*>(user_dataB);
  auto fn       = nb::cast<std::function<CVLsPrecSetupStdFnB>>(fn_table->lsprecsetupfnB);

  auto result = fn(t, y, yB, fyB, jokB, gammaB, nullptr);

  *jcurPtrB = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int cvode_lsprecsolvefnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsPrecSolveFnB>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::lsprecsolvefnB,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjactimessetupfnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesSetupFnB>, cvodea_user_supplied_fn_table,
    1>(&cvodea_user_supplied_fn_table::lsjactimessetupfnB,
       std::forward<Args>(args)...);
}

template<typename... Args>
inline int cvode_lsjactimesvecfnB_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<CVLsJacTimesVecFnB>, cvodea_user_supplied_fn_table,
    2>(&cvodea_user_supplied_fn_table::lsjactimesvecfnB,
       std::forward<Args>(args)...);
}

using CVLsLinSysStdFnB = std::tuple<int, sunbooleantype>(sunrealtype t, N_Vector y, 
  N_Vector yB, N_Vector fyB, SUNMatrix AB, sunbooleantype jokB, 
  sunrealtype gammaB, void *user_dataB, 
  N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

inline int cvode_lslinsysfnB_wrapper(sunrealtype t, N_Vector y, 
  N_Vector yB, N_Vector fyB, SUNMatrix AB, sunbooleantype jokB, 
  sunbooleantype *jcurB, sunrealtype gammaB, void *user_dataB, 
  N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  auto fn_table = static_cast<cvodea_user_supplied_fn_table*>(user_dataB);
  auto fn       = nb::cast<std::function<CVLsLinSysStdFnB>>(fn_table->lslinsysfnB);

  auto result = fn(t, y, yB, fyB, AB, jokB, gammaB, nullptr, tmp1B, tmp2B, tmp3B);

  *jcurB = std::get<1>(result);

  return std::get<0>(result);
}

//
// TODO(CJB): we can enable these functions with sundials v8.0.0
//            we need to add a int Ns argument to the callback like CVSensRhsFn has
//

// template<typename... Args>
// inline int cvode_fQBS_wrapper(Args... args)
// {
// }

// template<typename... Args>
// inline int cvode_lsjacfnBS_wrapper(Args... args)
// {
// }

// template<typename... Args>
// inline int cvode_lsprecsetupfnBS_wrapper(Args... args)
// {
// }

// template<typename... Args>
// inline int cvode_lsprecsolvefnBS_wrapper(Args... args)
// {
// }

// template<typename... Args>
// inline int cvode_lsjactimessetupfnBS_wrapper(Args... args)
// {
// }

// template<typename... Args>
// inline int cvode_lsjactimesvecfnBS_wrapper(Args... args)
// {
// }

// template<typename... Args>
// inline int cvode_lslinsysfnBS_wrapper(Args... args)
// {
// }

#endif
