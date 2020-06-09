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
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
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
  py_enum_<Spectral::Basis::Quadrature>

      // Wrapper for ExtentsAndTensorVolumeData
      py::class_<ExtentsAndTensorVolumeData>(m, "ExtentsAndTensorVolumeData")
          .def(py::init<std::vector<size_t>, std::vector<TensorComponent>>(),
               py::arg("extents"), py::arg("components"))
          .def("get_extents",
               [](const ExtentsAndTensorVolumeData& etvd) {
                 return etvd.extents;
               })
          .def("set_extents",
               [](ExtentsAndTensorVolumeData& etvd,
                  const std::vector<size_t>& v) { etvd.extents = v; })
          .def("set_tensor_components",
               [](ExtentsAndTensorVolumeData& etvd,
                  const std::vector<TensorComponent>& v_tc) {
                 etvd.tensor_components = v_tc;
               })
          .def("get_tensor_components",
               [](const ExtentsAndTensorVolumeData& etvd) {
                 return etvd.tensor_components;
               });

  py::class_<ElementVolumeData, ExtentsAndTensorVolumeData>(m,
                                                            "ElementVolumeData")
      .def(py::init<std::vector<size_t>, std::vector<TensorComponent>,
                    std::vector<Spectral::Basis>,
                    std::vector<Spectral::Quadrature>>(),
           py::arg("extents"), py::arg("components"), py::arg("basis"),
           py::arg("quadrature"))
      .def("get_extents",
           [](const ElementVolumeData& evd) { return etvd.extents; })
      .def("set_extents",
           [](ElementVolumeData& evd, const std::vector<size_t>& v) {
             etvd.extents = v;
           })
      .def(
          "set_tensor_components",
          [](ElementVolumeData& evd, const std::vector<TensorComponent>& v_tc) {
            etvd.tensor_components = v_tc;
          })
      .def("get_tensor_components",
           [](const ElementVolumeData& evd) { return evd.tensor_components; })
      .def("get_basis", [](const ElementVolumeData& evd) { return evd.basis })
      .def("get_quadrature",
           [](const ElementVolumeData& evd) { return evd.quadrature });
}

}  // namespace py_bindings
