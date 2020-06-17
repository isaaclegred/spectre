// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>
#include <string>
#include <utility>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TensorData.hpp"
#include "PythonBindings/BoundChecks.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/StdHelpers.hpp"

namespace py = pybind11;

namespace py_bindings {
void bind_tensordata(py::module& m) {  // NOLINT
  // Wrapper for TensorComponent
  py::class_<TensorComponent>(m, "TensorComponent")
      .def(py::init<std::string, DataVector>(), py::arg("name"),
           py::arg("data"))
      .def_readwrite("name", &TensorComponent::name)
      .def_readwrite("data", &TensorComponent::data)
      .def("__str__", [](const TensorComponent& tc) { return get_output(tc); })
      .def("__repr__",
           [](const TensorComponent& tc) { return get_output(tc); });

  // Wrapper for ExtentsAndTensorVolumeData
  py::class_<ExtentsAndTensorVolumeData>(m, "ExtentsAndTensorVolumeData")
      .def(py::init<std::vector<size_t>, std::vector<TensorComponent>>(),
           py::arg("extents"), py::arg("components"))
      .def_readwrite("extents", &ExtentsAndTensorVolumeData::extents)
      .def_readwrite("tensor_components",
                     &ExtentsAndTensorVolumeData::tensor_components)
      .def("__str__",
           [](const ExtentsAndTensorVolumeData& etvd) {
             std::string etvd_str = "(" + get_output(etvd.extents) + "," +
                                    get_output(etvd.tensor_components) + ")";
             return etvd_str;
           })
      .def("__repr__", [](const ExtentsAndTensorVolumeData& etvd) {
        std::string etvd_str = "(" + get_output(etvd.extents) + "," +
                               get_output(etvd.tensor_components) + ")";
        return etvd_str;
      });
}
}  // namespace py_bindings
