// #ifndef _NVECTOR_MANY_VECTOR_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_N_VectorContent_ManyVector =
  nb::class_<_N_VectorContent_ManyVector>(m, "_N_VectorContent_ManyVector", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "N_VNew_ManyVector",
  [](sunindextype num_subvectors, std::vector<N_Vector> vec_array_1d,
     SUNContext sunctx) -> N_Vector
  {
    auto N_VNew_ManyVector_adapt_arr_ptr_to_std_vector =
      [](sunindextype num_subvectors, std::vector<N_Vector> vec_array_1d,
         SUNContext sunctx) -> N_Vector
    {
      N_Vector* vec_array_1d_ptr = reinterpret_cast<N_Vector*>(
        vec_array_1d.empty() ? nullptr : vec_array_1d.data());

      auto lambda_result = N_VNew_ManyVector(num_subvectors, vec_array_1d_ptr,
                                             sunctx);
      return lambda_result;
    };

    return N_VNew_ManyVector_adapt_arr_ptr_to_std_vector(num_subvectors,
                                                         vec_array_1d, sunctx);
  },
  nb::arg("num_subvectors"), nb::arg("vec_array_1d"), nb::arg("sunctx"),
  nb::rv_policy::reference);

m.def("N_VGetSubvector_ManyVector", N_VGetSubvector_ManyVector, nb::arg("v"),
      nb::arg("vec_num"), nb::rv_policy::reference);

m.def("N_VGetNumSubvectors_ManyVector", N_VGetNumSubvectors_ManyVector,
      nb::arg("v"));

m.def("N_VGetSubvectorLocalLength_ManyVector",
      N_VGetSubvectorLocalLength_ManyVector, nb::arg("v"), nb::arg("vec_num"));

m.def("N_VEnableFusedOps_ManyVector", N_VEnableFusedOps_ManyVector,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableLinearCombination_ManyVector",
      N_VEnableLinearCombination_ManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleAddMulti_ManyVector", N_VEnableScaleAddMulti_ManyVector,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableDotProdMulti_ManyVector", N_VEnableDotProdMulti_ManyVector,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableLinearSumVectorArray_ManyVector",
      N_VEnableLinearSumVectorArray_ManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleVectorArray_ManyVector",
      N_VEnableScaleVectorArray_ManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableConstVectorArray_ManyVector",
      N_VEnableConstVectorArray_ManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormVectorArray_ManyVector",
      N_VEnableWrmsNormVectorArray_ManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormMaskVectorArray_ManyVector",
      N_VEnableWrmsNormMaskVectorArray_ManyVector, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableDotProdMultiLocal_ManyVector",
      N_VEnableDotProdMultiLocal_ManyVector, nb::arg("v"), nb::arg("tf"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
