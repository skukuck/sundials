// #ifndef _SUNDIALS_CONTEXT_H
//
// #ifdef __cplusplus
// #endif
//

m.def(
  "SUNContext_Create",
  [](SUNComm comm) -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNContext>>>
  {
    auto SUNContext_Create_adapt_modifiable_immutable_to_return =
      [](SUNComm comm) -> std::tuple<SUNErrCode, std::shared_ptr<std::remove_pointer_t<SUNContext>>>
    {
      SUNContext sunctx_out_adapt_modifiable;

      SUNErrCode r = SUNContext_Create(comm, &sunctx_out_adapt_modifiable);
      return std::make_tuple(r, our_make_shared<std::remove_pointer_t<SUNContext>, sundials::SUNContextDeleter>(sunctx_out_adapt_modifiable));
    };

    return SUNContext_Create_adapt_modifiable_immutable_to_return(comm);
  },
  nb::arg("comm"));

m.def("SUNContext_GetLastError", SUNContext_GetLastError, nb::arg("sunctx"));

m.def("SUNContext_PeekLastError", SUNContext_PeekLastError, nb::arg("sunctx"));

m.def("SUNContext_PushErrHandler", SUNContext_PushErrHandler, nb::arg("sunctx"),
      nb::arg("err_fn"), nb::arg("err_user_data"));

m.def("SUNContext_PopErrHandler", SUNContext_PopErrHandler, nb::arg("sunctx"));

m.def("SUNContext_ClearErrHandlers", SUNContext_ClearErrHandlers,
      nb::arg("sunctx"));

m.def(
  "SUNContext_GetProfiler",
  [](SUNContext sunctx) -> std::tuple<SUNErrCode, SUNProfiler>
  {
    auto SUNContext_GetProfiler_adapt_modifiable_immutable_to_return =
      [](SUNContext sunctx) -> std::tuple<SUNErrCode, SUNProfiler>
    {
      SUNProfiler profiler_adapt_modifiable;

      SUNErrCode r = SUNContext_GetProfiler(sunctx, &profiler_adapt_modifiable);
      return std::make_tuple(r, profiler_adapt_modifiable);
    };

    return SUNContext_GetProfiler_adapt_modifiable_immutable_to_return(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::rv_policy::reference,
  nb::keep_alive<0, 1>());

m.def("SUNContext_SetProfiler", SUNContext_SetProfiler, nb::arg("sunctx"),
      nb::arg("profiler"));

m.def(
  "SUNContext_GetLogger",
  [](SUNContext sunctx) -> std::tuple<SUNErrCode, SUNLogger>
  {
    auto SUNContext_GetLogger_adapt_modifiable_immutable_to_return =
      [](SUNContext sunctx) -> std::tuple<SUNErrCode, SUNLogger>
    {
      SUNLogger logger_adapt_modifiable;

      SUNErrCode r = SUNContext_GetLogger(sunctx, &logger_adapt_modifiable);
      return std::make_tuple(r, logger_adapt_modifiable);
    };

    return SUNContext_GetLogger_adapt_modifiable_immutable_to_return(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::rv_policy::reference,
  nb::keep_alive<0, 1>());

m.def("SUNContext_SetLogger", SUNContext_SetLogger, nb::arg("sunctx"),
      nb::arg("logger"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
