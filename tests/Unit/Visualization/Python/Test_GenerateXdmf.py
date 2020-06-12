# Distributed under the MIT License.
# See LICENSE.txt for details.

from spectre.Visualization.GenerateXdmf import generate_xdmf

import spectre.Informer as spectre_informer
import unittest
import os
from spectre.IO import H5 as spectre_h5
import spectre.DataStructures as ds


class TestGenerateXdmf(unittest.TestCase):
    def setUp(self):
        # We create a vol file using the python wrapper.
        self.file_name = os.path.join(spectre_informer.unit_test_path(),
                                      "Visualization/Python/VolTestData.h5")
        print(self.file_name)
        if os.path.isfile(self.file_name):
            os.remove(self.file_name)

        h5_file = spectre_h5.H5File(file_name=self.file_name,
                                    append_to_file=True)
        tensor_components_and_coords = [
            ds.DataVector([8.9, 7.6, 3.9, 2.1, 18.9, 17.6, 13.9, 12.1]),
            ds.DataVector([0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0]),
            ds.DataVector([0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0]),
            ds.DataVector([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]),
            ds.DataVector([-78.9, -7.6, -1.9, 1.9, 6.3, 1.7, 9.8, 1.2]),
            ds.DataVector([-7.9, 17.6, 1.9, -8.1, -6.3, 21.7, 6.8, -0.2]),
            ds.DataVector([1.0, 6.6, 28.0, -8.1, -26.3, 32.7, 26.8, -3.2])
        ]
        observation_ids = [8435087234, 6785087234]
        observation_values = [0.0, 1.0]
        grid_names = ["[[2,3,4]]", "[[5,6,7]]"]
        # Insert .vol file to h5 file
        h5_file.insert_vol("/element_data", 0)
        vol_file = h5_file.get_vol(path="/element_data")
        # Set TensorComponent and ExtentsAndTensorVolumeData to
        # be written
        basis = ds.Legendre
        quad = ds.Gauss
        element_vol_obs1 = [
            ds.ElementVolumeData([2, 2, 2], [
                ds.TensorComponent(grid_names[0] + "/S",
                                   tensor_components_and_coords[0]),
                ds.TensorComponent(grid_names[0] + "/InertialCoordinates_x",
                                   tensor_components_and_coords[1]),
                ds.TensorComponent(grid_names[0] + "/InertialCoordinates_y",
                                   tensor_components_and_coords[2]),
                ds.TensorComponent(grid_names[0] + "/InertialCoordinates_z",
                                   tensor_components_and_coords[3]),
                ds.TensorComponent(grid_names[0] + "/T_x",
                                   tensor_components_and_coords[4]),
                ds.TensorComponent(grid_names[0] + "/T_y",
                                   tensor_components_and_coords[5]),
                ds.TensorComponent(grid_names[0] + "/T_z",
                                   tensor_components_and_coords[6])
            ], [basis, basis, basis], [quad, quad, quad])
            for observation_value in observation_values
        ]

        element_vol_obs2 = [
            ds.ElementVolumeData([2, 2, 2], [
                ds.TensorComponent(grid_names[1] + "/S",
                                   tensor_components_and_coords[0]),
                ds.TensorComponent(grid_names[1] + "/x-coord",
                                   tensor_components_and_coords[2]),
                ds.TensorComponent(grid_names[1] + "/y-coord",
                                   tensor_components_and_coords[3]),
                ds.TensorComponent(grid_names[1] + "/z-coord",
                                   tensor_components_and_coords[1]),
                ds.TensorComponent(grid_names[1] + "/T_x",
                                   tensor_components_and_coords[6]),
                ds.TensorComponent(grid_names[1] + "/T_y",
                                   tensor_components_and_coords[4]),
                ds.TensorComponent(grid_names[1] + "/T_z",
                                   tensor_components_and_coords[5])
            ], [basis, basis, basis], [quad, quad, quad])
            for observation_value in observation_values
        ]
        # Write extents and tensor volume data to volfile
        vol_file.write_volume_data(observation_ids[0], observation_values[0],
                                   [element_vol_obs1[0], element_vol_obs2[0]])
        vol_file.write_volume_data(observation_ids[1], observation_values[1],
                                   [element_vol_obs1[1], element_vol_obs2[1]])

    def tearDown(self):
        if os.path.isfile(self.file_name):
            os.remove(self.file_name)

    def test_generate_xdmf(self):
        data_file_prefix = os.path.join(spectre_informer.unit_test_path(),
                                        'Visualization/Python', 'VolTestData')
        output_filename = 'Test_GenerateXdmf_output'
        if os.path.isfile(output_filename + '.xmf'):
            os.remove(output_filename + '.xmf')

        generate_xdmf(file_prefix=data_file_prefix,
                      output_filename=output_filename,
                      start_time=0.,
                      stop_time=1.,
                      stride=1,
                      coordinates='InertialCoordinates')

        # The script is quite opaque right now, so we only test that we can run
        # it and it produces output without raising an error. To test more
        # details, we should refactor the script into smaller units.
        self.assertTrue(os.path.isfile(output_filename + '.xmf'))
        os.remove(output_filename + '.xmf')


if __name__ == '__main__':
    unittest.main(verbosity=2)
