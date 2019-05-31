
// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <hdf5.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "IO/H5/Header.hpp"
#include "IO/H5/Object.hpp"
#include "IO/H5/OpenGroup.hpp"
#include "IO/H5/Python/Numpy.hpp"
#include "IO/H5/Python/ToNumpy.hpp"
#include "IO/H5/Version.hpp"
#include "IO/H5/VolumeData.hpp"
#include "PythonBindings/VectorPyList.hpp"

namespace bp = boost::python;

namespace py_bindings {
void bind_h5vol() {
  // Wrapper for basic H5VolumeData operations
  bp::class_<h5::VolumeData, boost::noncopyable>("H5Vol", bp::no_init)
      .def("extension", &h5::VolumeData::extension)
      // The method extension() is static
      .staticmethod("extension")

      .def("get_header", +[](h5::VolumeData& V) { return V.get_header(); })

      .def("get_version", +[](h5::VolumeData& V) { return V.get_version(); })

      .def("list_observation_ids",
           +[](h5::VolumeData& V) {
             return std_vector_to_py_list<size_t>(V.list_observation_ids());
           })
      .def("get_observation_value",
           +[](h5::VolumeData& V, size_t observation_id) {
             return V.get_observation_value(observation_id);
           })
      .def("list_grids",
           +[](h5::VolumeData& V, size_t observation_id) {
             return std_vector_to_py_list<std::string>(
                 V.list_grids(observation_id));
           })
      .def("list_tensor_components",
           +[](h5::VolumeData& V, size_t observation_id,
               const std::string& grid_name) {
             return std_vector_to_py_list<std::string>(
                 V.list_tensor_components(observation_id, grid_name));
           })
      .def("get_tensor_component",
           +[](h5::VolumeData& V, size_t observation_id,
               const std::string& grid_name,
               const std::string& tensor_component) {
             return V.get_tensor_component(observation_id, grid_name,
                                           tensor_component);
           })
      .def("get_extents", +[](h5::VolumeData& V, size_t observation_id,
                              const std::string& grid_name) {
        return std_vector_to_py_list(V.get_extents(observation_id, grid_name));
      });
}
}  // namespace py_bindings
