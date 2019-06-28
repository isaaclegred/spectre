# Distributed under the MIT License.
# See LICENSE.txt for details.

from spectre import DataStructures
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
        self.file_name_w = "pythontest.h5"
        # The other tests require external data, this file will not be changed
        self.file_name_r = Informer.unit_test_path() + "IO/VolTestData.h5"
        if os.path.isfile(self.file_name_w):
            os.remove(self.file_name_w)

    def tearDown(self):
        if os.path.isfile(self.file_name_w):
            os.remove(self.file_name_w)

    # Testing the VolumeData Insert Function
    def test_insert_vol(self):
        h5_file = spectre_h5.H5File(self.file_name_w, 1)
        h5_file.insert_vol("/element_data", 0)
        vol_file = h5_file.get_vol("/element_data")
        self.assertEqual(vol_file.get_version(), 0)
        h5_file.close()

    # Test the header was generated correctly
    def test_vol_get_header(self):
        h5_file = spectre_h5.H5File(self.file_name_w, 1)
        h5_file.insert_vol("/element_data", 0)
        vol_file = h5_file.get_vol("/element_data")
        self.assertEqual(vol_file.get_header()[0:20], "#\n# File created on ")
        h5_file.close()

    # The other tests require external data, so they use the included h5 file
    # Test the observation ids and values are correctly retrived
    def test_observation_id(self):
        h5_file = spectre_h5.H5File(self.file_name_r, 1)
        vol_file = h5_file.get_vol("/element_data")
        obs_ids = vol_file.list_observation_ids()
        expected_obs_ids = [836167400470363061, 1288279752782410135]
        expected_obs_values = [0.03, 0.02]
        self.assertItemsEqual(obs_ids, expected_obs_ids)
        for i, obs_id in enumerate(expected_obs_ids):
            self.assertEqual(vol_file.get_observation_value(obs_id),
                             expected_obs_values[i])
        h5_file.close()

    # Test to make sure information about the computation elements was found
    def test_grids(self):
        h5_file = spectre_h5.H5File(self.file_name_r, 1)
        vol_file = h5_file.get_vol("/element_data")
        obs_id = vol_file.list_observation_ids()[0]
        grid_names =  vol_file.get_grid_names(obs_id)
        expected_grid_names  = ['[B0,(L0I0,L0I0,L0I0)]']
        self.assertEqual(grid_names, expected_grid_names)
        extents = vol_file.get_extents(obs_id)
        expected_extents = [[2, 2, 2]]
        self.assertEqual(extents, expected_extents)
        h5_file.close()

    # Test that the tensor components, and tensor data  are retrieved correctly
    def test_tensor_compoenents(self):
        h5_file = spectre_h5.H5File(self.file_name_r, 1)
        vol_file = h5_file.get_vol("/element_data")
        obs_id = vol_file.list_observation_ids()[0]
        tensor_comps = vol_file.list_tensor_components(obs_id)
        expected_tensor_comps = ['Phi_x', 'Phi_y', 'Phi_z', 'Error(Phi)_x',
                                 'Error(Phi)_y', 'Error(Phi)_z',
                                 'InertialCoordinates_x',
                                 'InertialCoordinates_y',
                                 'InertialCoordinates_z']
        self.assertItemsEqual(tensor_comps, expected_tensor_comps)
        tensor_data = vol_file.get_tensor_component(
            obs_id, 'Phi_x')
        expected_tensor_data = DataStructures.DataVector([1.0, 1.0, 1.0, 1.0,
                                                        1.0, 1.0, 1.0, 1.0
                                                        ])
        self.assertEqual(tensor_data, expected_tensor_data)
        h5_file.close()


if __name__ == '__main__':
    unittest.main(verbosity=2)
