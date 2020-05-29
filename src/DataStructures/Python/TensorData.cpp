// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "DataStructures/Tensor/TensorData.hpp"
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>
#include <string>
#include <utility>

#include "DataStructures/DataVector.hpp"
#include "PythonBindings/BoundChecks.hpp"
#include "Utilities/GetOutput.hpp"

namespace py = pybind11;

namespace py_bindings {
void bind_tensordata(py::module& m) {  // NOLINT
  // Wrapper for TensorComponent
  py::class_<TensorComponent>(m, "TensorComponent")
      .def(py::init<std::string, DataVector>(), py::arg("name"),
           py::arg("data"))
      .def("set_name",
           [](TensorComponent& tc, const std::string& n) { tc.name = n; })
      .def("get_name", [](const TensorComponent& tc) { return tc.name; })
      .def("set_data",
           [](TensorComponent& tc, const DataVector& d) { tc.data = d; })
      .def("get_data", [](const TensorComponent& tc) { return tc.data; })
      .def("__str__", [](const TensorComponent& tc) { return get_output(tc); })
      .def("__repr__",
           [](const TensorComponent& tc) { return get_output(tc); });

  // Wrapper for ExtentsAndTensorVolumeData
  py::class_<ExtentsAndTensorVolumeData>(m, "ExtentsAndTensorVolumeData")
      .def(py::init<std::vector<size_t>, std::vector<TensorComponent>>(),
           py::arg("extents"), py::arg("components"))
      .def("get_extents",
           [](const ExtentsAndTensorVolumeData& etvd) { return etvd.extents; })
      .def("set_extents",
           [](ExtentsAndTensorVolumeData& etvd, const std::vector<size_t>& v) {
             etvd.extents = v;
           })
      .def("set_tensor_components",
           [](ExtentsAndTensorVolumeData& etvd,
              const std::vector<TensorComponent>& v_tc) {
             etvd.tensor_components = v_tc;
           })
      .def("get_tensor_components", [](const ExtentsAndTensorVolumeData& etvd) {
        return etvd.tensor_components;
      });
}
}  // namespace py_bindings
