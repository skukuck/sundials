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

#include "sundials4py.hpp"

#include <sundials/sundials_core.hpp>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes.hpp>
#include <cvodes/cvodes_ls.h>
#include <cvodes/cvodes_proj.h>

#include "cvodes/cvodes_impl.h"
#include "cvodes_usersupplied.hpp"

#include "sundials_adjointcheckpointscheme_impl.h"

namespace sundials4py {

using namespace sundials::experimental;

#define BIND_CVODE_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)       \
  m.def(                                                               \
    #NAME,                                                             \
    [](void* cv_mem, std::function<std::remove_pointer_t<FN_TYPE>> fn) \
    {                                                                  \
      auto fntable    = get_cvode_fn_table(cv_mem);                    \
      fntable->MEMBER = nb::cast(fn);                                  \
      if (fn) { return NAME(cv_mem, &WRAPPER); }                       \
      else { return NAME(cv_mem, nullptr); }                           \
    },                                                                 \
    __VA_ARGS__)

#define BIND_CVODE_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2, \
                             MEMBER2, WRAPPER2, ...)                      \
  m.def(                                                                  \
    #NAME,                                                                \
    [](void* cv_mem, std::function<std::remove_pointer_t<FN_TYPE1>> fn1,  \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                \
    {                                                                     \
      auto fntable     = get_cvode_fn_table(cv_mem);                      \
      fntable->MEMBER1 = nb::cast(fn1);                                   \
      fntable->MEMBER2 = nb::cast(fn2);                                   \
      if (fn1) { return NAME(cv_mem, WRAPPER1, WRAPPER2); }               \
      else { return NAME(cv_mem, nullptr, WRAPPER2); }                    \
    },                                                                    \
    __VA_ARGS__)

#define BIND_CVODEB_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)                 \
  m.def(                                                                          \
    #NAME,                                                                        \
    [](void* cv_mem, int which, std::function<std::remove_pointer_t<FN_TYPE>> fn) \
    {                                                                             \
      void* user_data = nullptr;                                                  \
      CVodeGetUserDataB(cv_mem, which, &user_data);                               \
      if (!user_data)                                                             \
        throw sundials4py::error_returned(                                        \
          "Failed to get Python function table from CVODE memory");               \
      auto fntable    = static_cast<cvodea_user_supplied_fn_table*>(user_data);   \
      fntable->MEMBER = nb::cast(fn);                                             \
      if (fn) { return NAME(cv_mem, which, &WRAPPER); }                           \
      else { return NAME(cv_mem, which, nullptr); }                               \
    },                                                                            \
    __VA_ARGS__)

#define BIND_CVODEB_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,       \
                              MEMBER2, WRAPPER2, ...)                            \
  m.def(                                                                         \
    #NAME,                                                                       \
    [](void* cv_mem, int which,                                                  \
       std::function<std::remove_pointer_t<FN_TYPE1>> fn1,                       \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                       \
    {                                                                            \
      void* user_data = nullptr;                                                 \
      CVodeGetUserDataB(cv_mem, which, &user_data);                              \
      if (!user_data)                                                            \
        throw sundials4py::error_returned(                                       \
          "Failed to get Python function table from CVODE memory");              \
      auto fntable     = static_cast<cvodea_user_supplied_fn_table*>(user_data); \
      fntable->MEMBER1 = nb::cast(fn1);                                          \
      fntable->MEMBER2 = nb::cast(fn2);                                          \
      if (fn1) { return NAME(cv_mem, which, WRAPPER1, WRAPPER2); }               \
      else { return NAME(cv_mem, which, nullptr, WRAPPER2); }                    \
    },                                                                           \
    __VA_ARGS__)

void bind_cvodes(nb::module_& m)
{
#include "cvodes_generated.hpp"

  nb::class_<CVodeView>(m, "CVodeView")
    .def("get", nb::overload_cast<>(&CVodeView::get, nb::const_),
         nb::rv_policy::reference);

  m.def(
    "CVodeCreate",
    [](int lmm, SUNContext sunctx)
    { return std::make_shared<CVodeView>(CVodeCreate(lmm, sunctx)); },
    nb::arg("lmm"), nb::arg("sunctx"), nb::keep_alive<0, 2>());

  m.def("CVodeInit",
        [](void* cv_mem, std::function<std::remove_pointer_t<CVRhsFn>> rhs,
           sunrealtype t0, N_Vector y0)
        {
          int cv_status = CVodeInit(cv_mem, cvode_f_wrapper, t0, y0);

          // Create the user-supplied function table to store the Python user functions
          auto cb_fns = cvode_user_supplied_fn_table_alloc();

          static_cast<CVodeMem>(cv_mem)->python = cb_fns;

          // Smuggle the user-supplied function table into callback wrappers through the user_data pointer
          cv_status = CVodeSetUserData(cv_mem, cv_mem);
          if (cv_status != CV_SUCCESS)
          {
            free(cb_fns);
            throw sundials4py::error_returned(
              "Failed to set user data in CVODE memory");
          }

          // Finally, set the RHS function
          cb_fns->f = nb::cast(rhs);

          return cv_status;
        });

  m.def("CVodeRootInit",
        [](void* cv_mem, int nrtfn,
           std::function<std::remove_pointer_t<CVRootFn>> fn)
        {
          auto fntable    = get_cvode_fn_table(cv_mem);
          fntable->rootfn = nb::cast(fn);
          return CVodeRootInit(cv_mem, nrtfn, &cvode_rootfn_wrapper);
        });

  m.def("CVodeQuadInit",
        [](void* cv_mem, std::function<std::remove_pointer_t<CVQuadRhsFn>> fQ,
           N_Vector yQ0)
        {
          auto fntable = get_cvode_fn_table(cv_mem);
          fntable->fQ  = nb::cast(fQ);
          return CVodeQuadInit(cv_mem, &cvode_fQ_wrapper, yQ0);
        });

  BIND_CVODE_CALLBACK(CVodeWFtolerances, CVEwtFn, ewtn, cvode_ewtfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("efun").none());

  BIND_CVODE_CALLBACK(CVodeSetNlsRhsFn, CVRhsFn, fNLS, cvode_nlsrhsfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("f").none());

  BIND_CVODE_CALLBACK(CVodeSetJacFn, CVLsJacFn, lsjacfn, cvode_lsjacfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("jac").none());

  BIND_CVODE_CALLBACK2(CVodeSetPreconditioner, CVLsPrecSetupFn, lsprecsetupfn,
                       cvode_lsprecsetupfn_wrapper, CVLsPrecSolveFn,
                       lsprecsolvefn, cvode_lsprecsolvefn_wrapper,
                       nb::arg("cvode_mem"), nb::arg("pset").none(),
                       nb::arg("psolve").none());

  BIND_CVODE_CALLBACK2(CVodeSetJacTimes, CVLsJacTimesSetupFn, lsjactimessetupfn,
                       cvode_lsjactimessetupfn_wrapper, CVLsJacTimesVecFn,
                       lsjactimesvecfn, cvode_lsjactimesvecfn_wrapper,
                       nb::arg("cvode_mem"), nb::arg("jtsetup").none(),
                       nb::arg("jtimes").none());

  BIND_CVODE_CALLBACK(CVodeSetLinSysFn, CVLsLinSysFn, lslinsysfn,
                      cvode_lslinsysfn_wrapper, nb::arg("cvode_mem"),
                      nb::arg("linsys").none());

  BIND_CVODE_CALLBACK(CVodeSetJacTimesRhsFn, CVLsJacTimesVecFn, lsjacrhsfn,
                      cvode_lsjacrhsfn_wrapper, nb::arg("cvode_mem"),
                      nb::arg("jtimesRhsFn").none());

  BIND_CVODE_CALLBACK(CVodeSetProjFn, CVProjFn, projfn, cvode_projfn_wrapper,
                      nb::arg("cvode_mem"), nb::arg("pfun").none());

  m.def("CVodeQuadSensInit",
        [](void* cv_mem, std::function<CVQuadSensRhsStdFn> fQS,
           std::vector<N_Vector> yQS0)
        {
          auto fntable = get_cvode_fn_table(cv_mem);
          fntable->fQS = nb::cast(fQS);
          return CVodeQuadSensInit(cv_mem, cvode_fQS_wrapper, yQS0.data());
        });

  m.def("CVodeSensInit",
        [](void* cv_mem, int Ns, int ism, std::function<CVSensRhsStdFn> fS,
           std::vector<N_Vector> yS0)
        {
          auto fntable = get_cvode_fn_table(cv_mem);
          fntable->fS  = nb::cast(fS);
          return CVodeSensInit(cv_mem, Ns, ism, cvode_fS_wrapper, yS0.data());
        });

  m.def("CVodeSensInit1",
        [](void* cv_mem, int Ns, int ism,
           std::function<std::remove_pointer_t<CVSensRhs1Fn>> fS1,
           std::vector<N_Vector> yS0)
        {
          auto fntable = get_cvode_fn_table(cv_mem);
          fntable->fS1 = nb::cast(fS1);
          return CVodeSensInit1(cv_mem, Ns, ism, cvode_fS1_wrapper, yS0.data());
        });

  ///
  // CVODES Adjoint Bindings
  ///

  m.def("CVodeInitB",
        [](void* cv_mem, int which,
           std::function<std::remove_pointer_t<CVRhsFnB>> fB, sunrealtype tB0,
           N_Vector yB0)
        {
          int cv_status = CVodeInitB(cv_mem, which, cvode_fB_wrapper, tB0, yB0);

          // Create the user-supplied function table to store the Python user functions
          auto cb_fns = cvodea_user_supplied_fn_table_alloc();

          // Store the function table in the python member of CVodeMem (if needed for global access)
          static_cast<CVodeMem>(cv_mem)->python = cb_fns;

          // Set user_data to cv_mem for compatibility (if needed)
          cv_status = CVodeSetUserDataB(cv_mem, which, cv_mem);
          if (cv_status != CV_SUCCESS)
          {
            free(cb_fns);
            throw sundials4py::error_returned(
              "Failed to set user data in CVODE memory");
          }

          // Finally, set the RHS function
          cb_fns->fB = nb::cast(fB);

          return cv_status;
        });

  m.def("CVodeQuadInitB",
        [](void* cv_mem, int which,
           std::function<std::remove_pointer_t<CVQuadRhsFnB>> fQB, N_Vector yQBO)
        {
          auto fntable = get_cvodea_fn_table(cv_mem);
          fntable->fQB = nb::cast(fQB);
          return CVodeQuadInitB(cv_mem, which, cvode_fQB_wrapper, yQBO);
        });

  BIND_CVODEB_CALLBACK(CVodeSetJacFnB, CVLsJacFnB, lsjacfnB,
                       cvode_lsjacfnB_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("jacB").none());

  BIND_CVODEB_CALLBACK2(CVodeSetPreconditionerB, CVLsPrecSetupFnB, lsprecsetupfnB,
                        cvode_lsprecsetupfnB_wrapper, CVLsPrecSolveFnB,
                        lsprecsolvefnB, cvode_lsprecsolvefnB_wrapper,
                        nb::arg("cv_mem"), nb::arg("which"),
                        nb::arg("psetB").none(), nb::arg("psolveB").none());

  BIND_CVODEB_CALLBACK2(CVodeSetJacTimesB, CVLsJacTimesSetupFnB,
                        lsjactimessetupfnB, cvode_lsjactimessetupfnB_wrapper,
                        CVLsJacTimesVecFnB, lsjactimesvecfnB,
                        cvode_lsjactimesvecfnB_wrapper, nb::arg("cv_mem"),
                        nb::arg("which"), nb::arg("jsetupB").none(),
                        nb::arg("jtimesB").none());

  BIND_CVODEB_CALLBACK(CVodeSetLinSysFnB, CVLsLinSysFnB, lslinsysfnB,
                       cvode_lslinsysfnB_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("linsysB").none());

  m.def("CVodeQuadInitBS",
        [](void* cv_mem, int which, std::function<CVQuadRhsStdFnBS> fQBS,
           N_Vector yQBO)
        {
          auto fntable  = get_cvodea_fn_table(cv_mem);
          fntable->fQBS = nb::cast(fQBS);
          return CVodeQuadInitBS(cv_mem, which, cvode_fQBS_wrapper, yQBO);
        });

  BIND_CVODEB_CALLBACK(CVodeSetJacFnBS, CVLsJacStdFnBS, lsjacfnBS,
                       cvode_lsjacfnBS_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("jacBS").none());

  BIND_CVODEB_CALLBACK2(CVodeSetPreconditionerBS, CVLsPrecSetupStdFnBS,
                        lsprecsetupfnBS, cvode_lsprecsetupfnBS_wrapper,
                        CVLsPrecSolveStdFnBS, lsprecsolvefnBS,
                        cvode_lsprecsolvefnBS_wrapper, nb::arg("cv_mem"),
                        nb::arg("which"), nb::arg("psetBS").none(),
                        nb::arg("psolveBS").none());

  BIND_CVODEB_CALLBACK2(CVodeSetJacTimesBS, CVLsJacTimesSetupStdFnBS,
                        lsjactimessetupfnBS, cvode_lsjactimessetupfnBS_wrapper,
                        CVLsJacTimesVecStdFnBS, lsjactimesvecfnBS,
                        cvode_lsjactimesvecfnBS_wrapper, nb::arg("cv_mem"),
                        nb::arg("which"), nb::arg("jsetupBS").none(),
                        nb::arg("jtimesBS").none());

  BIND_CVODEB_CALLBACK(CVodeSetLinSysFnBS, CVLsLinSysStdFnBS, lslinsysfnBS,
                       cvode_lslinsysfnBS_wrapper, nb::arg("cv_mem"),
                       nb::arg("which"), nb::arg("linsysBS").none());
}

} // namespace sundials4py
