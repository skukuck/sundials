#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_core.h>

#include "sundials4py_types.hpp"

namespace nb = nanobind;

namespace sundials4py {

void bind_nvector_serial(nb::module_& m)
{
#include "nvector_serial_generated.hpp"
}

} // namespace sundials4py
