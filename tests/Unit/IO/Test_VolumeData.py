# Distributed under the MIT License.
# See LICENSE.txt for details.

from spectre import DataStructures as ds
import spectre.IO.H5 as spectre_h5
from spectre import Informer
import unittest
import numpy as np
import os
import numpy.testing as npt


class TestIOH5VolumeData(unittest.TestCase):
    # Test Fixtures
    def setUp(self):
        # The First tests involve inserting vol files, the h5 file for
        # these will be deleted and recreated for each test as needed.
        self.file_name_w = os.path.join(Informer.unit_test_path(),
                                        "IO/pythontest.h5")
        # The other tests use a volume data file written using
        # the write_volume_data() function
        self.file_name_r = os.path.join(Informer.unit_test_path(),
                                        "IO/VolTestData.h5")
        if os.path.isfile(self.file_name_w):
            os.remove(self.file_name_w)
        if os.path.isfile(self.file_name_r):
            os.remove(self.file_name_r)

        h5_file = spectre_h5.H5File(file_name=self.file_name_r,
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
        observation_values = [7.0, 1.3]
        grid_names = ["[[2,3,4]]", "[[5,6,7]]"]

        # Insert .vol file to h5 file
        h5_file.insert_vol("/element_data", 0)
        vol_file = h5_file.get_vol(path="/element_data")
        # Set TensorComponent and ExtentsAndTensorVolumeData to
        # be written
        extents_and_tensor_vol_obs1 = [
            ds.ExtentsAndTensorVolumeData([2, 2, 2], [
                ds.TensorComponent(grid_names[0] + "/S",
                                   tensor_components_and_coords[0]),
                ds.TensorComponent(grid_names[0] + "/x-coord",
                                   tensor_components_and_coords[1]),
                ds.TensorComponent(grid_names[0] + "/y-coord",
                                   tensor_components_and_coords[2]),
                ds.TensorComponent(grid_names[0] + "/z-coord",
                                   tensor_components_and_coords[3]),
                ds.TensorComponent(grid_names[0] + "/T_x",
                                   tensor_components_and_coords[4]),
                ds.TensorComponent(grid_names[0] + "/T_y",
                                   tensor_components_and_coords[5]),
                ds.TensorComponent(grid_names[0] + "/T_z",
                                   tensor_components_and_coords[6])
            ]) for observation_value in observation_values
        ]

        extents_and_tensor_vol_obs2 = [
            ds.ExtentsAndTensorVolumeData([2, 2, 2], [
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
            ]) for observation_value in observation_values
        ]
        # Write extents and tensor volume data to volfile
        vol_file.write_volume_data(
            observation_ids[0], observation_values[0],
            [extents_and_tensor_vol_obs1[0], extents_and_tensor_vol_obs2[0]])
        vol_file.write_volume_data(
            observation_ids[1], observation_values[1],
            [extents_and_tensor_vol_obs1[1], extents_and_tensor_vol_obs2[1]])
        self.h5_file_r = h5_file

    def tearDown(self):
        if os.path.isfile(self.file_name_w):
            os.remove(self.file_name_w)
        if os.path.isfile(self.file_name_r):
            os.remove(self.file_name_r)

    # Testing the VolumeData Insert Function
    def test_insert_vol(self):
        h5_file = spectre_h5.H5File(file_name=self.file_name_w,
                                    append_to_file=True)
        h5_file.insert_vol(path="/element_data", version=0)
        vol_file = h5_file.get_vol(path="/element_data")
        self.assertEqual(vol_file.get_version(), 0)
        h5_file.close()

    # Test the header was generated correctly
    def test_vol_get_header(self):
        h5_file = spectre_h5.H5File(file_name=self.file_name_w,
                                    append_to_file=True)
        h5_file.insert_vol(path="/element_data", version=0)
        vol_file = h5_file.get_vol(path="/element_data")
        self.assertEqual(vol_file.get_header()[0:20], "#\n# File created on ")
        h5_file.close()

    # Using volume data H5 file to do rest of the tests
    # Test that observation ids and values are retrieved correctly
    def test_observation_id(self):
        h5_file = self.h5_file_r
        vol_file = h5_file.get_vol(path="/element_data")
        obs_ids = set(vol_file.list_observation_ids())
        expected_obs_ids = set([8435087234, 6785087234])
        expected_obs_values = {8435087234: 7.0, 6785087234: 1.3}
        self.assertEqual(obs_ids, expected_obs_ids)
        for obs_id in expected_obs_ids:
            self.assertEqual(
                vol_file.get_observation_value(observation_id=obs_id),
                expected_obs_values[obs_id])
        h5_file.close()

    # Test to make sure information about the computation elements was found
    def test_grids(self):
        h5_file = self.h5_file_r
        vol_file = h5_file.get_vol("/element_data")
        obs_id = vol_file.list_observation_ids()[0]
        grid_names = vol_file.get_grid_names(observation_id=obs_id)
        for grid_name in grid_names:
            str(grid_name)
        expected_grid_names = ["[[2,3,4]]", "[[5,6,7]]"]
        self.assertEqual(grid_names, expected_grid_names)
        extents = vol_file.get_extents(observation_id=obs_id)
        expected_extents = [[2, 2, 2], [2, 2, 2]]
        self.assertEqual(extents, expected_extents)
        h5_file.close()

    # Test that the tensor components, and tensor data  are retrieved correctly
    def test_tensor_components(self):
        h5_file = self.h5_file_r
        vol_file = h5_file.get_vol(path="/element_data")
        obs_id = vol_file.list_observation_ids()[0]
        tensor_comps = set(
            vol_file.list_tensor_components(observation_id=obs_id))
        expected_tensor_comps = set(
            ['S', 'x-coord', 'y-coord', 'z-coord', 'T_x', 'T_y', 'T_z'])
        self.assertEqual(tensor_comps, expected_tensor_comps)
        expected_S_tensor_data = np.array(
            [8.9, 7.6, 3.9, 2.1, 18.9, 17.6, 13.9, 12.1])
        expected_xcoord_tensor_data = np.array(
            [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0])
        expected_ycoord_tensor_data = np.array(
            [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0])
        expected_zcoord_tensor_data = np.array(
            [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
        expected_Tx_tensor_data = np.array(
            [-78.9, -7.6, -1.9, 1.9, 6.3, 1.7, 9.8, 1.2])
        expected_Ty_tensor_data = np.array(
            [-7.9, 17.6, 1.9, -8.1, -6.3, 21.7, 6.8, -0.2])
        expected_Tz_tensor_data = np.array(
            [1.0, 6.6, 28.0, -8.1, -26.3, 32.7, 26.8, -3.2])
        # Checking whether two numpy arrays are "almost equal" is easy, so
        # we convert everything to numpy arrays for comparison.
        S_tensor_data = np.asarray(
            vol_file.get_tensor_component(observation_id=obs_id,
                                          tensor_component='S'))
        npt.assert_array_almost_equal(S_tensor_data[0:8],
                                      expected_S_tensor_data)
        xcoord_tensor_data = np.asarray(
            vol_file.get_tensor_component(observation_id=obs_id,
                                          tensor_component='x-coord'))
        npt.assert_array_almost_equal(xcoord_tensor_data[0:8],
                                      expected_xcoord_tensor_data)
        ycoord_tensor_data = np.asarray(
            vol_file.get_tensor_component(observation_id=obs_id,
                                          tensor_component='y-coord'))
        npt.assert_array_almost_equal(ycoord_tensor_data[0:8],
                                      expected_ycoord_tensor_data)
        zcoord_tensor_data = np.asarray(
            vol_file.get_tensor_component(observation_id=obs_id,
                                          tensor_component='z-coord'))
        npt.assert_array_almost_equal(zcoord_tensor_data[0:8],
                                      expected_zcoord_tensor_data)
        Tx_tensor_data = np.asarray(
            vol_file.get_tensor_component(observation_id=obs_id,
                                          tensor_component='T_x'))
        npt.assert_array_almost_equal(Tx_tensor_data[0:8],
                                      expected_Tx_tensor_data)
        Ty_tensor_data = np.asarray(
            vol_file.get_tensor_component(observation_id=obs_id,
                                          tensor_component='T_y'))
        npt.assert_array_almost_equal(Ty_tensor_data[0:8],
                                      expected_Ty_tensor_data)
        Tz_tensor_data = np.asarray(
            vol_file.get_tensor_component(observation_id=obs_id,
                                          tensor_component='T_z'))
        npt.assert_array_almost_equal(Tz_tensor_data[0:8],
                                      expected_Tz_tensor_data)
        h5_file.close()


if __name__ == '__main__':
    unittest.main(verbosity=2)
