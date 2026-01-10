// #ifndef _NVECTOR_CUDA_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClass_N_VectorContent_Cuda =
  nb::class_<_N_VectorContent_Cuda>(m, "_N_VectorContent_Cuda", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "N_VNewEmpty_Cuda",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VNewEmpty_Cuda_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result = N_VNewEmpty_Cuda(sunctx);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VNewEmpty_Cuda_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());

m.def(
  "N_VNew_Cuda",
  [](sunindextype length,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VNew_Cuda_adapt_return_type_to_shared_ptr =
      [](sunindextype length,
         SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result = N_VNew_Cuda(length, sunctx);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VNew_Cuda_adapt_return_type_to_shared_ptr(length, sunctx);
  },
  nb::arg("length"), nb::arg("sunctx"), "nb::keep_alive<0, 2>()",
  nb::keep_alive<0, 2>());

m.def(
  "N_VNewManaged_Cuda",
  [](sunindextype length,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VNewManaged_Cuda_adapt_return_type_to_shared_ptr =
      [](sunindextype length,
         SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result = N_VNewManaged_Cuda(length, sunctx);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VNewManaged_Cuda_adapt_return_type_to_shared_ptr(length, sunctx);
  },
  nb::arg("length"), nb::arg("sunctx"), "nb::keep_alive<0, 2>()",
  nb::keep_alive<0, 2>());

m.def(
  "N_VNewWithMemHelp_Cuda",
  [](sunindextype length, sunbooleantype use_managed_mem, SUNMemoryHelper helper,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VNewWithMemHelp_Cuda_adapt_return_type_to_shared_ptr =
      [](sunindextype length, sunbooleantype use_managed_mem,
         SUNMemoryHelper helper,
         SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result = N_VNewWithMemHelp_Cuda(length, use_managed_mem,
                                                  helper, sunctx);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VNewWithMemHelp_Cuda_adapt_return_type_to_shared_ptr(length,
                                                                  use_managed_mem,
                                                                  helper, sunctx);
  },
  nb::arg("length"), nb::arg("use_managed_mem"), nb::arg("helper"),
  nb::arg("sunctx"), "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());

m.def(
  "N_VMake_Cuda",
  [](sunindextype length, sundials4py::Array1d h_vdata_1d,
     sundials4py::Array1d d_vdata_1d,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VMake_Cuda_adapt_arr_ptr_to_std_vector =
      [](sunindextype length, sundials4py::Array1d h_vdata_1d,
         sundials4py::Array1d d_vdata_1d, SUNContext sunctx) -> N_Vector
    {
      sunrealtype* h_vdata_1d_ptr = h_vdata_1d.size() == 0 ? nullptr
                                                           : h_vdata_1d.data();
      sunrealtype* d_vdata_1d_ptr = d_vdata_1d.size() == 0 ? nullptr
                                                           : d_vdata_1d.data();

      auto lambda_result = N_VMake_Cuda(length, h_vdata_1d_ptr, d_vdata_1d_ptr,
                                        sunctx);
      return lambda_result;
    };
    auto N_VMake_Cuda_adapt_return_type_to_shared_ptr =
      [&N_VMake_Cuda_adapt_arr_ptr_to_std_vector](sunindextype length,
                                                  sundials4py::Array1d h_vdata_1d,
                                                  sundials4py::Array1d d_vdata_1d,
                                                  SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result =
        N_VMake_Cuda_adapt_arr_ptr_to_std_vector(length, h_vdata_1d, d_vdata_1d,
                                                 sunctx);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VMake_Cuda_adapt_return_type_to_shared_ptr(length, h_vdata_1d,
                                                        d_vdata_1d, sunctx);
  },
  nb::arg("length"), nb::arg("h_vdata_1d"), nb::arg("d_vdata_1d"),
  nb::arg("sunctx"), "nb::keep_alive<0, 4>()", nb::keep_alive<0, 4>());

m.def(
  "N_VMakeManaged_Cuda",
  [](sunindextype length, sundials4py::Array1d vdata_1d,
     SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
  {
    auto N_VMakeManaged_Cuda_adapt_arr_ptr_to_std_vector =
      [](sunindextype length, sundials4py::Array1d vdata_1d,
         SUNContext sunctx) -> N_Vector
    {
      sunrealtype* vdata_1d_ptr = vdata_1d.size() == 0 ? nullptr
                                                       : vdata_1d.data();

      auto lambda_result = N_VMakeManaged_Cuda(length, vdata_1d_ptr, sunctx);
      return lambda_result;
    };
    auto N_VMakeManaged_Cuda_adapt_return_type_to_shared_ptr =
      [&N_VMakeManaged_Cuda_adapt_arr_ptr_to_std_vector](sunindextype length,
                                                         sundials4py::Array1d vdata_1d,
                                                         SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      auto lambda_result =
        N_VMakeManaged_Cuda_adapt_arr_ptr_to_std_vector(length, vdata_1d, sunctx);

      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        lambda_result);
    };

    return N_VMakeManaged_Cuda_adapt_return_type_to_shared_ptr(length, vdata_1d,
                                                               sunctx);
  },
  nb::arg("length"), nb::arg("vdata_1d"), nb::arg("sunctx"),
  "nb::keep_alive<0, 3>()", nb::keep_alive<0, 3>());

m.def(
  "N_VSetHostArrayPointer_Cuda",
  [](sundials4py::Array1d h_vdata_1d, N_Vector v)
  {
    auto N_VSetHostArrayPointer_Cuda_adapt_arr_ptr_to_std_vector =
      [](sundials4py::Array1d h_vdata_1d, N_Vector v)
    {
      sunrealtype* h_vdata_1d_ptr = h_vdata_1d.size() == 0 ? nullptr
                                                           : h_vdata_1d.data();

      N_VSetHostArrayPointer_Cuda(h_vdata_1d_ptr, v);
    };

    N_VSetHostArrayPointer_Cuda_adapt_arr_ptr_to_std_vector(h_vdata_1d, v);
  },
  nb::arg("h_vdata_1d"), nb::arg("v"));

m.def(
  "N_VSetDeviceArrayPointer_Cuda",
  [](sundials4py::Array1d d_vdata_1d, N_Vector v)
  {
    auto N_VSetDeviceArrayPointer_Cuda_adapt_arr_ptr_to_std_vector =
      [](sundials4py::Array1d d_vdata_1d, N_Vector v)
    {
      sunrealtype* d_vdata_1d_ptr = d_vdata_1d.size() == 0 ? nullptr
                                                           : d_vdata_1d.data();

      N_VSetDeviceArrayPointer_Cuda(d_vdata_1d_ptr, v);
    };

    N_VSetDeviceArrayPointer_Cuda_adapt_arr_ptr_to_std_vector(d_vdata_1d, v);
  },
  nb::arg("d_vdata_1d"), nb::arg("v"));

m.def("N_VIsManagedMemory_Cuda", N_VIsManagedMemory_Cuda, nb::arg("x"));

m.def("N_VSetKernelExecPolicy_Cuda", N_VSetKernelExecPolicy_Cuda, nb::arg("x"),
      nb::arg("stream_exec_policy"), nb::arg("reduce_exec_policy"));

m.def("N_VCopyToDevice_Cuda", N_VCopyToDevice_Cuda, nb::arg("v"));

m.def("N_VCopyFromDevice_Cuda", N_VCopyFromDevice_Cuda, nb::arg("v"));

m.def("N_VGetHostArrayPointer_Cuda", N_VGetHostArrayPointer_Cuda, nb::arg("x"));

m.def("N_VEnableFusedOps_Cuda", N_VEnableFusedOps_Cuda, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableLinearCombination_Cuda", N_VEnableLinearCombination_Cuda,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleAddMulti_Cuda", N_VEnableScaleAddMulti_Cuda, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableDotProdMulti_Cuda", N_VEnableDotProdMulti_Cuda, nb::arg("v"),
      nb::arg("tf"));

m.def("N_VEnableLinearSumVectorArray_Cuda", N_VEnableLinearSumVectorArray_Cuda,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleVectorArray_Cuda", N_VEnableScaleVectorArray_Cuda,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableConstVectorArray_Cuda", N_VEnableConstVectorArray_Cuda,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormVectorArray_Cuda", N_VEnableWrmsNormVectorArray_Cuda,
      nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableWrmsNormMaskVectorArray_Cuda",
      N_VEnableWrmsNormMaskVectorArray_Cuda, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableScaleAddMultiVectorArray_Cuda",
      N_VEnableScaleAddMultiVectorArray_Cuda, nb::arg("v"), nb::arg("tf"));

m.def("N_VEnableLinearCombinationVectorArray_Cuda",
      N_VEnableLinearCombinationVectorArray_Cuda, nb::arg("v"), nb::arg("tf"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
