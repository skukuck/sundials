#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>

#include <nvector/nvector_manyvector.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_nvector.hpp>

#include "sundials4py_types.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_nvector_manyvector(nb::module_& m)
{
#include "nvector_manyvector_generated.hpp"
}

} // namespace sundials4py
