// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "tests/Unit/TestingFramework.hpp"

#include <algorithm>
#include <boost/iterator/transform_iterator.hpp>
#include <hdf5.h>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TensorData.hpp"
#include "ErrorHandling/Assert.hpp"
#include "ErrorHandling/Error.hpp"
#include "IO/Connectivity.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/Header.hpp"
#include "IO/H5/Helpers.hpp"
#include "IO/H5/Version.hpp"
#include "Utilities/Numeric.hpp"

namespace h5 {
void check_volume_data(
    const VolumeData& volume_file, const size_t temporal_id,
    const std::vector<ExtentsAndTensorVolumeData>& expected_data) {
  // Check that the temporal_id is in the file
  auto read_observation_ids = volume_file.list_observation_ids();
  CHECK(std::find(read_observation_ids.begin(), read_observation_ids.end(),
                  temporal_id) != read_observation_ids.end());

  // Check the dimension of the grids is correct
  size_t expected_dim = expected_data.front().extents.size();
  CHECK(expected_dim == volume_file.get_dimension());

  // Check that all of the tensor components were written as expected
  auto read_component_names = volume_file.list_tensor_components(temporal_id);
  const auto get_component_name = [](const auto& component) noexcept {
    ASSERT(component.name.find_last_of('/') != std::string::npos,
           "The expected format of the tensor component names is "
           "'GROUP_NAME/COMPONENT_NAME' but could not find a '/' in '"
               << component.name << "'.");
    return component.name.substr(component.name.find_last_of('/') + 1);
  };

  const std::vector<std::string> expected_component_names(
      boost::make_transform_iterator(
          expected_data.front().tensor_components.begin(), get_component_name),
      boost::make_transform_iterator(
          expected_data.front().tensor_components.end(), get_component_name));

  // If the two lists are the same size, and one contains another, they are
  // equal
  CHECK(read_component_names.size() == expected_component_names.size());
  REQUIRE(alg::all_of(read_component_names,
                      [&expected_component_names](const std::string& name) {
                        return alg::found(expected_component_names, name);
                      }));
  // Check that each element's data was written correctly
  auto expected_element_names = volume_file.get_grid_names(temporal_id);
  CHECK(expected_element_names.size() == expected_data.size());
  // We need to know how many points in the data correspond to the elements
  // written before a particular element, say E, so we can determine where in
  // the contiguous DataVectors to look for E's data.
  // `read_points_by_element[i]` stores the number of points which were written
  // "by element `i`" so that the data for element `i` is found starting at
  // `read_points_by_element[i]` in all datavectors.

  // Helper Function to get number of points on a particular grid
  auto accumulate_extents = [](std::vector<size_t> grid_extents) {
    return alg::accumulate(grid_extents, 1, std::multiplies<>{});
  };

  auto read_extents = volume_file.get_extents(temporal_id);
  std::vector<size_t> element_num_points(
      boost::make_transform_iterator(read_extents.begin(), accumulate_extents),
      boost::make_transform_iterator(read_extents.end(), accumulate_extents));
  auto read_points_by_element = [&element_num_points]() {
    std::vector<size_t> read_points(element_num_points.size());
    read_points[0] = 0;
    size_t index = 1;
    while (index != element_num_points.size()) {
      read_points[index] =
          read_points[index - 1] + element_num_points[index - 1];
      index++;
    }
    return read_points;
  }();

  for (auto element_data : expected_data) {
    const auto& first_tensor_name = element_data.tensor_components.front().name;
    ASSERT(first_tensor_name.find_last_of('/') != std::string::npos,
           "The expected format of the tensor component names is "
           "'GROUP_NAME/COMPONENT_NAME' but could not find a '/' in '"
               << first_tensor_name << "'.");
    const auto element_name =
        first_tensor_name.substr(0, first_tensor_name.find_last_of('/'));
    auto position = std::find(expected_element_names.begin(),
                              expected_element_names.end(), element_name);
    // Make sure that the component was found
    CHECK(position != expected_element_names.end());
    size_t index = static_cast<size_t>(
        std::distance(expected_element_names.begin(), position));
    // Check the extents of the Element
    CHECK(read_extents[index] == element_data.extents);
    // Check the TensorData on the element
    for (auto&& component : element_data.tensor_components) {
      auto read_component_data = volume_file.get_tensor_component(
          temporal_id, get_component_name(component));
      CHECK(DataVector(&read_component_data[read_points_by_element[index]],
                       element_num_points[index]) == component.data);
    }
  }
}
}  // namespace h5
