// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <boost/python.hpp>

namespace py_bindings {
void bind_infofrombuild();
}  // namespace py_bindings

BOOST_PYTHON_MODULE(_Informer) {
  Py_Initialize();
  py_bindings::bind_infofrombuild();
}
