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

#include <idas/idas.h>
#include <idas/idas.hpp>
#include <idas/idas_ls.h>

#include "idas/idas_impl.h"

#include "idas_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

#define BIND_IDA_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)          \
  m.def(                                                                \
    #NAME,                                                              \
    [](void* ida_mem, std::function<std::remove_pointer_t<FN_TYPE>> fn) \
    {                                                                   \
      auto fn_table    = get_idas_fn_table(ida_mem);                    \
      fn_table->MEMBER = nb::cast(fn);                                  \
      if (fn) { return NAME(ida_mem, &WRAPPER); }                       \
      else { return NAME(ida_mem, nullptr); }                           \
    },                                                                  \
    __VA_ARGS__)

#define BIND_IDA_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2,   \
                           MEMBER2, WRAPPER2, ...)                        \
  m.def(                                                                  \
    #NAME,                                                                \
    [](void* ida_mem, std::function<std::remove_pointer_t<FN_TYPE1>> fn1, \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)                \
    {                                                                     \
      auto fn_table     = get_idas_fn_table(ida_mem);                     \
      fn_table->MEMBER1 = nb::cast(fn1);                                  \
      fn_table->MEMBER2 = nb::cast(fn2);                                  \
      if (fn1) { return NAME(ida_mem, WRAPPER1, WRAPPER2); }              \
      else { return NAME(ida_mem, nullptr, WRAPPER2); }                   \
    },                                                                    \
    __VA_ARGS__)

#define BIND_IDAB_CALLBACK(NAME, FN_TYPE, MEMBER, WRAPPER, ...)                    \
  m.def(                                                                           \
    #NAME,                                                                         \
    [](void* ida_mem, int which, std::function<std::remove_pointer_t<FN_TYPE>> fn) \
    {                                                                              \
      void* user_data  = nullptr;                                                  \
      auto fn_table    = get_idasa_fn_table(ida_mem, which);                       \
      fn_table->MEMBER = nb::cast(fn);                                             \
      if (fn) { return NAME(ida_mem, which, &WRAPPER); }                           \
      else { return NAME(ida_mem, which, nullptr); }                               \
    },                                                                             \
    __VA_ARGS__)

#define BIND_IDAB_CALLBACK2(NAME, FN_TYPE1, MEMBER1, WRAPPER1, FN_TYPE2, \
                            MEMBER2, WRAPPER2, ...)                      \
  m.def(                                                                 \
    #NAME,                                                               \
    [](void* ida_mem, int which,                                         \
       std::function<std::remove_pointer_t<FN_TYPE1>> fn1,               \
       std::function<std::remove_pointer_t<FN_TYPE2>> fn2)               \
    {                                                                    \
      void* user_data   = nullptr;                                       \
      auto fn_table     = get_idasa_fn_table(ida_mem, which);            \
      fn_table->MEMBER1 = nb::cast(fn1);                                 \
      fn_table->MEMBER2 = nb::cast(fn2);                                 \
      if (fn1) { return NAME(ida_mem, which, WRAPPER1, WRAPPER2); }      \
      else { return NAME(ida_mem, which, nullptr, WRAPPER2); }           \
    },                                                                   \
    __VA_ARGS__)

namespace sundials4py {

void bind_idas(nb::module_& m)
{
#include "idas_generated.hpp"

  nb::class_<IDAView>(m, "IDAView")
    .def("get", nb::overload_cast<>(&IDAView::get, nb::const_),
         nb::rv_policy::reference);

  m.def(
    "IDACreate",
    [](SUNContext sunctx)
    { return std::make_shared<IDAView>(IDACreate(sunctx)); },
    nb::arg("sunctx"), nb::keep_alive<0, 1>());

  m.def("IDAInit",
        [](void* ida_mem, std::function<std::remove_pointer_t<IDAResFn>> res,
           sunrealtype t0, N_Vector yy0, N_Vector yp0)
        {
          int ida_status = IDAInit(ida_mem, idas_res_wrapper, t0, yy0, yp0);

          auto cb_fns = idas_user_supplied_fn_table_alloc();
          static_cast<IDAMem>(ida_mem)->python = cb_fns;

          ida_status = IDASetUserData(ida_mem, ida_mem);
          if (ida_status != IDA_SUCCESS)
          {
            free(cb_fns);
            throw sundials4py::error_returned(
              "Failed to set user data in IDAS memory");
          }

          cb_fns->res = nb::cast(res);

          return ida_status;
        });

  m.def("IDARootInit",
        [](void* ida_mem, int nrtfn,
           std::function<std::remove_pointer_t<IDARootStdFn>> fn)
        {
          auto fn_table    = get_idas_fn_table(ida_mem);
          fn_table->rootfn = nb::cast(fn);
          return IDARootInit(ida_mem, nrtfn, &idas_rootfn_wrapper);
        });

  m.def("IDAQuadInit",
        [](void* ida_mem,
           std::function<std::remove_pointer_t<IDAQuadRhsFn>> resQ, N_Vector yQ0)
        {
          auto fn_table  = get_idas_fn_table(ida_mem);
          fn_table->resQ = nb::cast(resQ);
          return IDAQuadInit(ida_mem, &idas_resQ_wrapper, yQ0);
        });

  BIND_IDA_CALLBACK(IDAWFtolerances, IDAEwtFn, ewtn, idas_ewtfn_wrapper,
                    nb::arg("ida_mem"), nb::arg("efun").none());

  BIND_IDA_CALLBACK(IDASetNlsResFn, IDAResFn, resNLS, idas_nlsresfn_wrapper,
                    nb::arg("ida_mem"), nb::arg("res").none());

  BIND_IDA_CALLBACK(IDASetJacFn, IDALsJacFn, lsjacfn, idas_lsjacfn_wrapper,
                    nb::arg("ida_mem"), nb::arg("jac").none());

  BIND_IDA_CALLBACK2(IDASetPreconditioner, IDALsPrecSetupFn, lsprecsetupfn,
                     idas_lsprecsetupfn_wrapper, IDALsPrecSolveFn, lsprecsolvefn,
                     idas_lsprecsolvefn_wrapper, nb::arg("ida_mem"),
                     nb::arg("pset").none(), nb::arg("psolve").none());

  BIND_IDA_CALLBACK2(IDASetJacTimes, IDALsJacTimesSetupFn, lsjactimessetupfn,
                     idas_lsjactimessetupfn_wrapper, IDALsJacTimesVecFn,
                     lsjactimesvecfn, idas_lsjactimesvecfn_wrapper,
                     nb::arg("ida_mem"), nb::arg("jtsetup").none(),
                     nb::arg("jtimes").none());

  BIND_IDA_CALLBACK(IDASetJacTimesResFn, IDALsJacTimesVecFn, lsjacresfn,
                    idas_lsjacresfn_wrapper, nb::arg("ida_mem"),
                    nb::arg("jtimesResFn").none());

  //
  // Sensitivity and quadrature sensitivity bindings
  //

  m.def("IDAQuadSensInit",
        [](void* ida_mem, std::function<IDAQuadSensRhsStdFn> resQS,
           std::vector<N_Vector> yQS0)
        {
          auto fn_table   = get_idas_fn_table(ida_mem);
          fn_table->resQS = nb::cast(resQS);
          return IDAQuadSensInit(ida_mem, idas_resQS_wrapper, yQS0.data());
        });

  m.def("IDASensInit",
        [](void* ida_mem, int Ns, int ism, std::function<IDASensResStdFn> resS,
           std::vector<N_Vector> yS0, std::vector<N_Vector> ypS0)
        {
          auto fn_table  = get_idas_fn_table(ida_mem);
          fn_table->resS = nb::cast(resS);
          return IDASensInit(ida_mem, Ns, ism, idas_resS_wrapper, yS0.data(),
                             ypS0.data());
        });

  ///
  // IDAS adjoint bindings
  ///

  m.def("IDAInitB",
        [](void* ida_mem, int which,
           std::function<std::remove_pointer_t<IDAResFnB>> resB,
           sunrealtype tB0, N_Vector yyB0, N_Vector ypB0)
        {
          int ida_status = IDAInitB(ida_mem, which, idas_resB_wrapper, tB0,
                                    yyB0, ypB0);

          auto cb_fns   = idasa_user_supplied_fn_table_alloc();
          auto idab_mem = static_cast<IDAMem>(IDAGetAdjIDABmem(ida_mem, which));
          idab_mem->python = cb_fns;

          ida_status = IDASetUserDataB(ida_mem, which, idab_mem);
          if (ida_status != IDA_SUCCESS)
          {
            free(cb_fns);
            throw sundials4py::error_returned(
              "Failed to set user data in IDAS memory");
          }

          cb_fns->resB = nb::cast(resB);
          return ida_status;
        });

  m.def("IDAQuadInitB",
        [](void* ida_mem, int which,
           std::function<std::remove_pointer_t<IDAQuadRhsFnB>> resQB,
           N_Vector yQBO)
        {
          auto fn_table   = get_idasa_fn_table(ida_mem, which);
          fn_table->resQB = nb::cast(resQB);
          return IDAQuadInitB(ida_mem, which, idas_resQB_wrapper, yQBO);
        });

  BIND_IDAB_CALLBACK(IDASetJacFnB, IDALsJacFnB, lsjacfnB, idas_lsjacfnB_wrapper,
                     nb::arg("ida_mem"), nb::arg("which"),
                     nb::arg("jacB").none());

  BIND_IDAB_CALLBACK2(IDASetPreconditionerB, IDALsPrecSetupFnB, lsprecsetupfnB,
                      idas_lsprecsetupfnB_wrapper, IDALsPrecSolveFnB,
                      lsprecsolvefnB, idas_lsprecsolvefnB_wrapper,
                      nb::arg("ida_mem"), nb::arg("which"),
                      nb::arg("psetB").none(), nb::arg("psolveB").none());

  BIND_IDAB_CALLBACK2(IDASetJacTimesB, IDALsJacTimesSetupFnB, lsjactimessetupfnB,
                      idas_lsjactimessetupfnB_wrapper, IDALsJacTimesVecFnB,
                      lsjactimesvecfnB, idas_lsjactimesvecfnB_wrapper,
                      nb::arg("ida_mem"), nb::arg("which"),
                      nb::arg("jsetupB").none(), nb::arg("jtimesB").none());

  m.def("IDAInitBS",
        [](void* ida_mem, int which, std::function<IDAResStdFnBS> resBS,
           sunrealtype tB0, N_Vector yyB0, N_Vector ypB0)
        {
          int ida_status = IDAInitBS(ida_mem, which, ida_resBS_wrapper, tB0,
                                     yyB0, ypB0);

          auto cb_fns   = idasa_user_supplied_fn_table_alloc();
          auto idab_mem = static_cast<IDAMem>(IDAGetAdjIDABmem(ida_mem, which));
          idab_mem->python = cb_fns;

          ida_status = IDASetUserDataB(ida_mem, which, idab_mem);
          if (ida_status != IDA_SUCCESS)
          {
            free(cb_fns);
            throw sundials4py::error_returned(
              "Failed to set user data in IDA memory");
          }

          cb_fns->resBS = nb::cast(resBS);
          return ida_status;
        });

  m.def("IDAQuadInitBS",
        [](void* ida_mem, int which, std::function<IDAQuadRhsStdFnBS> resQBS,
           N_Vector yQBO)
        {
          auto fn_table    = get_idasa_fn_table(ida_mem, which);
          fn_table->resQBS = nb::cast(resQBS);
          return IDAQuadInitBS(ida_mem, which, idas_resQBS_wrapper, yQBO);
        });

  BIND_IDAB_CALLBACK(IDASetJacFnBS, IDALsJacStdFnBS, lsjacfnBS,
                     idas_lsjacfnBS_wrapper, nb::arg("ida_mem"),
                     nb::arg("which"), nb::arg("jacBS").none());

  BIND_IDAB_CALLBACK2(IDASetPreconditionerBS, IDALsPrecSetupStdFnBS,
                      lsprecsetupfnBS, idas_lsprecsetupfnBS_wrapper,
                      IDALsPrecSolveStdFnBS, lsprecsolvefnBS,
                      idas_lsprecsolvefnBS_wrapper, nb::arg("ida_mem"),
                      nb::arg("which"), nb::arg("psetBS").none(),
                      nb::arg("psolveBS").none());

  BIND_IDAB_CALLBACK2(IDASetJacTimesBS, IDALsJacTimesSetupStdFnBS,
                      lsjactimessetupfnBS, idas_lsjactimessetupfnBS_wrapper,
                      IDALsJacTimesVecStdFnBS, lsjactimesvecfnBS,
                      idas_lsjactimesvecfnBS_wrapper, nb::arg("ida_mem"),
                      nb::arg("which"), nb::arg("jsetupBS").none(),
                      nb::arg("jtimesBS").none());
}

} // namespace sundials4py
