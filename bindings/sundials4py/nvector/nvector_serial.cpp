#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>

// #include <sundials/sundials_core.hpp>
#include <sundials/sundials_nvector.hpp>
#include <nvector/nvector_serial.h>

#include "sundials/sundials_classview.hpp"
#include "sundials4py_types.hpp"

namespace nb = nanobind;
using sundials::experimental::NVectorDeleter;
using sundials::experimental::our_make_shared;

namespace sundials4py {

void bind_nvector_serial(nb::module_& m)
{
#include "nvector_serial_generated.hpp"

  m.def(
    "N_VNew_Serial",
    [](sunindextype length,
       std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx)
    {
      return our_make_shared<std::remove_pointer_t<N_Vector>, NVectorDeleter>(N_VNew_Serial(length, sunctx.get()));
    },
    nb::arg("vec_length"), nb::arg("sunctx"), nb::keep_alive<0, 2>());
}

} // namespace sundials4py
