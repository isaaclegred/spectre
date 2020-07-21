// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <algorithm>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

#include "Utilities/Algorithm.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TensorData.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/StdHelpers.hpp"

namespace py = pybind11;

namespace {
  double pi = M_PI;
  std::vector<ExtentsAndTensorVolumeData> get_sphere_data(double R, size_t theta_pts,
                                                          size_t phi_pts){
    DataVector thetas(theta_pts);
    DataVector phis(phi_pts);
    DataVector theta_inds(theta_pts);
    DataVector phi_inds(phi_pts);
    std::iota(theta_inds.begin(), theta_inds.end(), 0);
    std::iota(phi_inds.begin(), phi_inds.end(), 0);
    size_t num_points = theta_pts * phi_pts;
    alg::transform(theta_inds, thetas.begin(),
                   [theta_pts](const size_t& pt){
                     return static_cast<double>(pt)*pi/static_cast<double>(theta_pts);}
);
    alg::transform(phi_inds, phis.begin(),
                   [phi_pts](const size_t& pt){
                     return static_cast<double>(pt)*2.0*pi/static_cast<double>(phi_pts);});
  

  DataVector x_coords(num_points);
  DataVector y_coords(num_points);
  DataVector z_coords(num_points);
  size_t point = 0;
  for (auto theta : thetas){
      for (auto phi : phis){
        z_coords[point] = R * cos(theta);
        x_coords[point] = R * sin(theta) * cos(phi);
        y_coords[point] = R * sin(theta) * sin(phi);
        point++;
      }
  }
  return std::vector<ExtentsAndTensorVolumeData>{{{theta_pts, phi_pts},{
     TensorComponent("InertialCoordinates_x", x_coords),
        TensorComponent("InertialCoordinates_y", y_coords),
       TensorComponent("InertialCoordinates_z", z_coords)}}};
  }
}


namespace py_bindings {
void bind_tensordata(py::module& m) {  // NOLINT
  // Wrapper for TensorComponent
  py::class_<TensorComponent>(m, "TensorComponent")
      .def(py::init<std::string, DataVector>(), py::arg("name"),
           py::arg("data"))
      .def_readwrite("name", &TensorComponent::name)
      .def_readwrite("data", &TensorComponent::data)
      .def("__str__", get_output<TensorComponent>)
      .def("__repr__", get_output<TensorComponent>)
      // NOLINTNEXTLINE(misc-redundant-expression)
      .def(py::self == py::self)
      // NOLINTNEXTLINE(misc-redundant-expression)
      .def(py::self != py::self);

  // Wrapper for ExtentsAndTensorVolumeData
  py::class_<ExtentsAndTensorVolumeData>(m, "ExtentsAndTensorVolumeData")
      .def(py::init<std::vector<size_t>, std::vector<TensorComponent>>(),
           py::arg("extents"), py::arg("components"))
      .def_readwrite("extents", &ExtentsAndTensorVolumeData::extents)
      .def_readwrite("tensor_components",
                     &ExtentsAndTensorVolumeData::tensor_components)
      .def("__str__",
           [](const ExtentsAndTensorVolumeData& extents_and_data) {
             return "(" + get_output(extents_and_data.extents) + "," +
                    get_output(extents_and_data.tensor_components) + ")";
           })
      .def("__repr__", [](const ExtentsAndTensorVolumeData& extents_and_data) {
        return "(" + get_output(extents_and_data.extents) + "," +
          get_output(extents_and_data.tensor_components) + ")";});
  m.def("get_sphere_data", &get_sphere_data, py::arg("R"), py::arg("theta_pts"),
        py::arg("phi_pts"));


}
}  // namespace py_bindings
