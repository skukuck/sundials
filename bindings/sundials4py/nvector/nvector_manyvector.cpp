#include "sundials4py.hpp"

#include <nvector/nvector_manyvector.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_nvector.hpp>

namespace nb = nanobind;
using namespace sundials::experimental;
using namespace sundials::experimental;

namespace sundials4py {

void bind_nvector_manyvector(nb::module_& m)
{
#include "nvector_manyvector_generated.hpp"
}

} // namespace sundials4py
