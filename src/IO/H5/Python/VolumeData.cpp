// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>

#include "Utilities/Algorithm.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TensorData.hpp"
#include "IO/H5/VolumeData.hpp"

namespace py = pybind11;

namespace {
  double pi = M_PI;
  h5::VolumeData get_sphere_data(double r, size_t theta_pts, size_t phi_pts){
    DataVector thetas(theta_pts);
    DataVector phis(phi_pts);
    std::iota(thetas.begin(), thetas.end());
    std::iota(phis.begin(), phis.end());
    size_t num_points = theta_pts * phi_pts;
    alg::transform(thetas,
                   [theta_pts](const size_t& pt){
                     return static_cast<double>(pt)*pi/static_cast<double>(theta_pts);},
                   thetas.begin());
    alg::transform(phis,
                   [phi_pts](const size_t& pt){
                     return static_cast<double>(pt)*2.0*pi/static_cast<double>(phi_pts);},
                   phis.begin());
  }

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
  auto cells =  vis::detail::compute_cells(Sphere<2>(theta_pts, phi_pts));
}

namespace py_bindings {
void bind_h5vol(py::module& m) {  // NOLINT
  // Wrapper for basic H5VolumeData operations
  py::class_<h5::VolumeData>(m, "H5Vol")
      .def_static("extension", &h5::VolumeData::extension)
      .def("get_header", &h5::VolumeData::get_header)
    .def("get_sphere", )
      .def("get_version", &h5::VolumeData::get_version)
      .def("write_volume_data", &h5::VolumeData::write_volume_data)
      .def("list_observation_ids", &h5::VolumeData::list_observation_ids)
      .def("get_observation_value", &h5::VolumeData::get_observation_value,
           py::arg("observation_id"))
      .def("get_grid_names", &h5::VolumeData::get_grid_names,
           py::arg("observation_id"))
      .def("list_tensor_components", &h5::VolumeData::list_tensor_components,
           py::arg("observation_id"))
      .def("get_tensor_component", &h5::VolumeData::get_tensor_component,
           py::arg("observation_id"), py::arg("tensor_component"))
      .def("get_extents", &h5::VolumeData::get_extents,
           py::arg("observation_id"));
  m.def("offset_and_length_for_grid", &h5::offset_and_length_for_grid,
        py::arg("grid_name"), py::arg("all_grid_names"),
        py::arg("all_extents"));
}
}  // namespace py_bindings
