#include "sundials4py.hpp"

#include <nvector/nvector_cuda.h>
#include <sundials/sundials_nvector.hpp>

#include "sundials/sundials_classview.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_nvector_cuda(nb::module_& m)
{
#include "nvector_cuda_generated.hpp"
}

} // namespace sundials4py
