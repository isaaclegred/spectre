// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "tests/Unit/TestingFramework.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TensorData.hpp"
#include "ErrorHandling/Error.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/VolumeData.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/Numeric.hpp"
#include "tests/Unit/IO/IOHelpers.hpp"


SPECTRE_TEST_CASE("Unit.IO.H5.VolumeData", "[Unit][IO][H5]") {
  const std::string h5_file_name("Unit.IO.H5.VolumeData.h5");
  const uint32_t version_number = 4;
  if (file_system::check_if_file_exists(h5_file_name)) {
    file_system::rm(h5_file_name, true);
  }
  h5::H5File<h5::AccessType::ReadWrite> my_file(h5_file_name);
  const std::vector<DataVector> tensor_components_and_coords{
      {8.9, 7.6, 3.9, 2.1, 18.9, 17.6, 13.9, 12.1},
      {0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0},
      {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
      {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
      {-78.9, -7.6, -1.9, 8.1, 6.3, 8.7, 9.8, 0.2},
      {-7.9, 7.6, 1.9, -8.1, -6.3, 2.7, 6.8, -0.2},
      {17.9, 27.6, 21.9, -28.1, -26.3, 32.7, 26.8, -30.2}};
  const std::vector<size_t> observation_ids{8435087234, size_t(-1)};
  const std::vector<double> observation_values{8.0, 2.3};
  const std::vector<std::string> grid_names{"[[2,3,4]]", "[[5,6,7]]"};
  auto make_volume_data = [&tensor_components_and_coords,
                           &grid_names](const double observation_value) {
    std::string first_grid = grid_names.front();
    std::string last_grid = grid_names.back();
    return std::vector<ExtentsAndTensorVolumeData>{
        {{2, 2, 2},
         {TensorComponent{first_grid + "/S",
                          observation_value * tensor_components_and_coords[0]},
          TensorComponent{first_grid + "/x-coord",
                          observation_value * tensor_components_and_coords[1]},
          TensorComponent{first_grid + "/y-coord",
                          observation_value * tensor_components_and_coords[2]},
          TensorComponent{first_grid + "/z-coord",
                          observation_value * tensor_components_and_coords[3]},
          TensorComponent{first_grid + "/T_x",
                          observation_value * tensor_components_and_coords[4]},
          TensorComponent{first_grid + "/T_y",
                          observation_value * tensor_components_and_coords[5]},
          TensorComponent{
              first_grid + "/T_z",
              observation_value * tensor_components_and_coords[6]}}},
        // Second Element Data
        {{2, 2, 2},
         {TensorComponent{last_grid + "/S",
                          observation_value * tensor_components_and_coords[1]},
          TensorComponent{last_grid + "/x-coord",
                          observation_value * tensor_components_and_coords[0]},
          TensorComponent{last_grid + "/y-coord",
                          observation_value * tensor_components_and_coords[5]},
          TensorComponent{last_grid + "/z-coord",
                          observation_value * tensor_components_and_coords[3]},
          TensorComponent{last_grid + "/T_x",
                          observation_value * tensor_components_and_coords[6]},
          TensorComponent{last_grid + "/T_y",
                          observation_value * tensor_components_and_coords[4]},
          TensorComponent{
              last_grid + "/T_z",
              observation_value * tensor_components_and_coords[2]}}}};
  };
  {
    auto& volume_file =
        my_file.insert<h5::VolumeData>("/element_data", version_number);
    const auto write_to_file = [
      &volume_file, &tensor_components_and_coords, &grid_names, make_volume_data
    ](const size_t observation_id, const double observation_value) noexcept {
      std::string first_grid = grid_names.front();
      std::string last_grid = grid_names.back();
      volume_file.write_volume_data(
          observation_id, observation_value,
          std::vector<ExtentsAndTensorVolumeData>{
              {{2, 2, 2},
               {TensorComponent{
                    first_grid + "/S",
                    observation_value * tensor_components_and_coords[0]},
                TensorComponent{
                    first_grid + "/x-coord",
                    observation_value * tensor_components_and_coords[1]},
                TensorComponent{
                    first_grid + "/y-coord",
                    observation_value * tensor_components_and_coords[2]},
                TensorComponent{
                    first_grid + "/z-coord",
                    observation_value * tensor_components_and_coords[3]},
                TensorComponent{
                    first_grid + "/T_x",
                    observation_value * tensor_components_and_coords[4]},
                TensorComponent{
                    first_grid + "/T_y",
                    observation_value * tensor_components_and_coords[5]},
                TensorComponent{
                    first_grid + "/T_z",
                    observation_value * tensor_components_and_coords[6]}}},
              // Second Element Data
              {{2, 2, 2},
               {TensorComponent{
                    last_grid + "/S",
                    observation_value * tensor_components_and_coords[1]},
                TensorComponent{
                    last_grid + "/x-coord",
                    observation_value * tensor_components_and_coords[0]},
                TensorComponent{
                    last_grid + "/y-coord",
                    observation_value * tensor_components_and_coords[5]},
                TensorComponent{
                    last_grid + "/z-coord",
                    observation_value * tensor_components_and_coords[3]},
                TensorComponent{
                    last_grid + "/T_x",
                    observation_value * tensor_components_and_coords[6]},
                TensorComponent{
                    last_grid + "/T_y",
                    observation_value * tensor_components_and_coords[4]},
                TensorComponent{
                    last_grid + "/T_z",
                    observation_value * tensor_components_and_coords[2]}}}});
    };

    for (size_t i = 0; i < observation_ids.size(); ++i) {
      write_to_file(observation_ids[i], observation_values[i]);
    }
  }
  // Open the read volume file and check that the observation id and values are
  // correct.

  const auto& volume_file =
      my_file.get<h5::VolumeData>("/element_data", version_number);
  const auto read_observation_ids = volume_file.list_observation_ids();
  // Check that all of the observation_ids in the file were correctly written
  // ids (If it is in the file it was written correctly)
  CHECK(alg::all_of(read_observation_ids, [&observation_ids](const size_t id) {
    return alg::found(observation_ids, id);
  }));
  // Check that the volume data is correct
  for (size_t i = 0; i < observation_ids.size(); ++i) {
    const auto& target_data =
      make_volume_data(volume_file.get_observation_value(observation_ids[i]));
    // Check observation_ids[i] exsits in file and has correct data
    check_volume_data(volume_file, observation_ids[i], target_data);
  }
  if (file_system::check_if_file_exists(h5_file_name)) {
    // file_system::rm(h5_file_name, true);
  }
}

// [[OutputRegex, The expected format of the tensor component names is
// 'GROUP_NAME/COMPONENT_NAME' but could not find a '/' in]]
[[noreturn]] SPECTRE_TEST_CASE("Unit.IO.H5.VolumeData.ComponentFormat0",
                               "[Unit][IO][H5]") {
  ASSERTION_TEST();
#ifdef SPECTRE_DEBUG
  const std::string h5_file_name("Unit.IO.H5.VolumeData.ComponentFormat.h5");
  const uint32_t version_number = 4;
  if (file_system::check_if_file_exists(h5_file_name)) {
    file_system::rm(h5_file_name, true);
  }
  h5::H5File<h5::AccessType::ReadWrite> my_file(h5_file_name);
  auto& volume_file =
      my_file.insert<h5::VolumeData>("/element_data", version_number);
  volume_file.write_volume_data(100, 10.0,
                                std::vector<ExtentsAndTensorVolumeData>{
                                    {{2}, {TensorComponent{"S", {1.0, 2.0}}}}});
  ERROR("Failed to trigger ASSERT in an assertion test");
#endif
  // clang-format off
}

// [[OutputRegex, The expected format of the tensor component names is
// 'GROUP_NAME/COMPONENT_NAME' but could not find a '/' in]]
[[noreturn]]
SPECTRE_TEST_CASE("Unit.IO.H5.VolumeData.ComponentFormat1",
                               "[Unit][IO][H5]") {
  // clang-format on
  ASSERTION_TEST();
#ifdef SPECTRE_DEBUG
  const std::string h5_file_name("Unit.IO.H5.VolumeData.ComponentFormat1.h5");
  const uint32_t version_number = 4;
  if (file_system::check_if_file_exists(h5_file_name)) {
    file_system::rm(h5_file_name, true);
  }
  h5::H5File<h5::AccessType::ReadWrite> my_file(h5_file_name);
  auto& volume_file =
      my_file.insert<h5::VolumeData>("/element_data", version_number);
  volume_file.write_volume_data(
      100, 10.0,
      std::vector<ExtentsAndTensorVolumeData>{
          ExtentsAndTensorVolumeData({2}, {TensorComponent{"A/S", {1.0, 2.0}},
                                           TensorComponent{"S", {1.0, 2.0}}})});
  ERROR("Failed to trigger ASSERT in an assertion test");
#endif
  // clang-format off
}


// [[OutputRegex, Trying to write tensor component 'S' which already exists
// in HDF5 file in group 'element_data.vol/ObservationId100']]
[[noreturn]] SPECTRE_TEST_CASE("Unit.IO.H5.VolumeData.WriteTwice",
                               "[Unit][IO][H5]") {
  // clang-format on
  ASSERTION_TEST();
#ifdef SPECTRE_DEBUG
  const std::string h5_file_name("Unit.IO.H5.VolumeData.WriteTwice.h5");
  const uint32_t version_number = 4;
  if (file_system::check_if_file_exists(h5_file_name)) {
    file_system::rm(h5_file_name, true);
  }
  h5::H5File<h5::AccessType::ReadWrite> my_file(h5_file_name);
  auto& volume_file =
      my_file.insert<h5::VolumeData>("/element_data", version_number);
  volume_file.write_volume_data(100, 10.0,
                                std::vector<ExtentsAndTensorVolumeData>{
                                    {{2},
                                     {TensorComponent{"A/S", {1.0, 2.0}},
                                      TensorComponent{"A/S", {1.0, 2.0}}}}});
  ERROR("Failed to trigger ASSERT in an assertion test");
#endif
  // clang-format off
}
