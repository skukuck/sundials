// #ifndef _ARKODE_SPRKSTEP_H
//
// #ifdef __cplusplus
// #endif
//

m.def("SPRKStepSetMethod", SPRKStepSetMethod, nb::arg("arkode_mem"),
      nb::arg("sprk_storage"));

m.def("SPRKStepSetMethodName", SPRKStepSetMethodName, nb::arg("arkode_mem"),
      nb::arg("method"));

m.def(
  "SPRKStepGetCurrentMethod",
  [](void* arkode_mem)
    -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>>
  {
    auto SPRKStepGetCurrentMethod_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, ARKodeSPRKTable>
    {
      ARKodeSPRKTable sprk_storage_adapt_modifiable;

      int r = SPRKStepGetCurrentMethod(arkode_mem,
                                       &sprk_storage_adapt_modifiable);
      return std::make_tuple(r, sprk_storage_adapt_modifiable);
    };
    auto SPRKStepGetCurrentMethod_adapt_return_type_to_shared_ptr =
      [&SPRKStepGetCurrentMethod_adapt_modifiable_immutable_to_return](
        void* arkode_mem)
      -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>>
    {
      auto lambda_result =
        SPRKStepGetCurrentMethod_adapt_modifiable_immutable_to_return(arkode_mem);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<ARKodeSPRKTable>,
                                             ARKodeSPRKTableDeleter>(
                               std::get<1>(lambda_result)));
    };

    return SPRKStepGetCurrentMethod_adapt_return_type_to_shared_ptr(arkode_mem);
  },
  nb::arg("arkode_mem"), nb::rv_policy::reference);
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
