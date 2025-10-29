// #ifndef _SUNDIALS_ADAPTCONTROLLER_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumSUNAdaptController_Type =
  nb::enum_<SUNAdaptController_Type>(m, "SUNAdaptController_Type",
                                     nb::is_arithmetic(), "")
    .value("SUN_ADAPTCONTROLLER_NONE", SUN_ADAPTCONTROLLER_NONE, "")
    .value("SUN_ADAPTCONTROLLER_H", SUN_ADAPTCONTROLLER_H, "")
    .value("SUN_ADAPTCONTROLLER_MRI_H_TOL", SUN_ADAPTCONTROLLER_MRI_H_TOL, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClass_generic_SUNAdaptController_Ops =
  nb::class_<_generic_SUNAdaptController_Ops>(m,
                                              "_generic_SUNAdaptController_Ops",
                                              "")
    .def(nb::init<>()) // implicit default constructor
  ;

auto pyClass_generic_SUNAdaptController =
  nb::class_<_generic_SUNAdaptController>(m, "_generic_SUNAdaptController", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def("SUNAdaptController_GetType", SUNAdaptController_GetType, nb::arg("C"));

m.def("SUNAdaptController_EstimateStep", SUNAdaptController_EstimateStep,
      nb::arg("C"), nb::arg("h"), nb::arg("p"), nb::arg("dsm"), nb::arg("hnew"));

m.def("SUNAdaptController_EstimateStepTol", SUNAdaptController_EstimateStepTol,
      nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("P"),
      nb::arg("DSM"), nb::arg("dsm"), nb::arg("Hnew"), nb::arg("tolfacnew"));

m.def("SUNAdaptController_Reset", SUNAdaptController_Reset, nb::arg("C"));

m.def("SUNAdaptController_SetDefaults", SUNAdaptController_SetDefaults,
      nb::arg("C"));

m.def("SUNAdaptController_Write", SUNAdaptController_Write, nb::arg("C"),
      nb::arg("fptr"));

m.def("SUNAdaptController_SetErrorBias", SUNAdaptController_SetErrorBias,
      nb::arg("C"), nb::arg("bias"));

m.def("SUNAdaptController_UpdateH", SUNAdaptController_UpdateH, nb::arg("C"),
      nb::arg("h"), nb::arg("dsm"));

m.def("SUNAdaptController_UpdateMRIHTol", SUNAdaptController_UpdateMRIHTol,
      nb::arg("C"), nb::arg("H"), nb::arg("tolfac"), nb::arg("DSM"),
      nb::arg("dsm"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
