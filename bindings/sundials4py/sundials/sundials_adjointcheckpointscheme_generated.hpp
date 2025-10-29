// #ifndef _SUNADJOINT_CHECKPOINTSCHEME_H
//
// #ifdef __cplusplus
// #endif
//

m.def("SUNAdjointCheckpointScheme_NeedsSaving",
      SUNAdjointCheckpointScheme_NeedsSaving, nb::arg("check_scheme"),
      nb::arg("step_num"), nb::arg("stage_num"), nb::arg("t"),
      nb::arg("yes_or_no"));

m.def("SUNAdjointCheckpointScheme_InsertVector",
      SUNAdjointCheckpointScheme_InsertVector, nb::arg("check_scheme"),
      nb::arg("step_num"), nb::arg("stage_num"), nb::arg("t"), nb::arg("state"));

m.def(
  "SUNAdjointCheckpointScheme_LoadVector",
  [](SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
     suncountertype stage_num, sunbooleantype peek,
     sunrealtype* tout) -> std::tuple<SUNErrCode, N_Vector>
  {
    auto SUNAdjointCheckpointScheme_LoadVector_adapt_modifiable_immutable_to_return =
      [](SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
         suncountertype stage_num, sunbooleantype peek,
         sunrealtype* tout) -> std::tuple<SUNErrCode, N_Vector>
    {
      N_Vector out_adapt_modifiable;

      SUNErrCode r =
        SUNAdjointCheckpointScheme_LoadVector(check_scheme, step_num, stage_num,
                                              peek, &out_adapt_modifiable, tout);
      return std::make_tuple(r, out_adapt_modifiable);
    };

    return SUNAdjointCheckpointScheme_LoadVector_adapt_modifiable_immutable_to_return(check_scheme,
                                                                                      step_num,
                                                                                      stage_num,
                                                                                      peek,
                                                                                      tout);
  },
  nb::arg("check_scheme"), nb::arg("step_num"), nb::arg("stage_num"),
  nb::arg("peek"), nb::arg("tout"), nb::rv_policy::reference);

m.def("SUNAdjointCheckpointScheme_EnableDense",
      SUNAdjointCheckpointScheme_EnableDense, nb::arg("check_scheme"),
      nb::arg("on_or_off"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
