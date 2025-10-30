// #ifndef _SUNADJOINTCHECKPOINTSCHEME_FIXED_H
//
// #ifdef __cplusplus
// #endif
//

m.def(
  "SUNAdjointCheckpointScheme_Create_Fixed",
  [](SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, suncountertype interval,
     suncountertype estimate, sunbooleantype keep, SUNContext sunctx)
    -> std::tuple<SUNErrCode,
                  std::shared_ptr<std::remove_pointer_t<SUNAdjointCheckpointScheme>>>
  {
    auto SUNAdjointCheckpointScheme_Create_Fixed_adapt_modifiable_immutable_to_return =
      [](SUNDataIOMode io_mode, SUNMemoryHelper mem_helper,
         suncountertype interval, suncountertype estimate, sunbooleantype keep,
         SUNContext sunctx) -> std::tuple<SUNErrCode, SUNAdjointCheckpointScheme>
    {
      SUNAdjointCheckpointScheme check_scheme_ptr_adapt_modifiable;

      SUNErrCode r =
        SUNAdjointCheckpointScheme_Create_Fixed(io_mode, mem_helper, interval,
                                                estimate, keep, sunctx,
                                                &check_scheme_ptr_adapt_modifiable);
      return std::make_tuple(r, check_scheme_ptr_adapt_modifiable);
    };
    auto SUNAdjointCheckpointScheme_Create_Fixed_adapt_return_type_to_shared_ptr =
      [&SUNAdjointCheckpointScheme_Create_Fixed_adapt_modifiable_immutable_to_return](SUNDataIOMode io_mode,
                                                                                      SUNMemoryHelper
                                                                                        mem_helper,
                                                                                      suncountertype
                                                                                        interval,
                                                                                      suncountertype
                                                                                        estimate,
                                                                                      sunbooleantype
                                                                                        keep,
                                                                                      SUNContext sunctx)
      -> std::tuple<SUNErrCode,
                    std::shared_ptr<std::remove_pointer_t<SUNAdjointCheckpointScheme>>>
    {
      auto lambda_result =
        SUNAdjointCheckpointScheme_Create_Fixed_adapt_modifiable_immutable_to_return(io_mode,
                                                                                     mem_helper,
                                                                                     interval,
                                                                                     estimate,
                                                                                     keep,
                                                                                     sunctx);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<
                               std::remove_pointer_t<SUNAdjointCheckpointScheme>,
                               SUNAdjointCheckpointSchemeDeleter>(
                               std::get<1>(lambda_result)));
    };

    return SUNAdjointCheckpointScheme_Create_Fixed_adapt_return_type_to_shared_ptr(io_mode,
                                                                                   mem_helper,
                                                                                   interval,
                                                                                   estimate,
                                                                                   keep,
                                                                                   sunctx);
  },
  nb::arg("io_mode"), nb::arg("mem_helper"), nb::arg("interval"),
  nb::arg("estimate"), nb::arg("keep"), nb::arg("sunctx"),
  "nb::call_policy<sundials4py::returns_references_to<6, 1>>()",
  nb::rv_policy::reference,
  nb::call_policy<sundials4py::returns_references_to<6, 1>>());
// #ifdef __cplusplus
//
// #endif
//
// #endif
