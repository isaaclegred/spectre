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

    # The other tests require external data
    def test_observation_id(self):
        h5_file = spectre_h5.H5File(self.file_name_r, 1)
        vol_file = h5_file.get_vol("/element_data")
        obs_id = vol_file.list_observation_ids()[0]
        self.assertEqual(obs_id, 16436106908031328247L)
        self.assertEqual(vol_file.get_observation_value(obs_id), .01)
        h5_file.close()

    def test_grids(self):
        h5_file = spectre_h5.H5File(self.file_name_r, 1)
        vol_file = h5_file.get_vol("/element_data")
        obs_id = vol_file.list_observation_ids()[0]
        grids = vol_file.list_grids(obs_id)
        target_grids = []
        for first in xrange(0, 2):
            for second in xrange(0, 2):
                for third in xrange(0, 2):
                    target_grids.append(
                        '[B0,(L1I'+str(first)+',L1I'+str(second)+',L1I' +
                        str(third)+')]')
        self.assertEqual(grids, target_grids)
        extents = vol_file.get_extents(obs_id, grids[0])
        self.assertEqual(extents, [2, 2, 2])
        h5_file.close()

    def test_tensor_compoenents(self):
        h5_file = spectre_h5.H5File(self.file_name_r, 1)
        vol_file = h5_file.get_vol("/element_data")
        obs_id = vol_file.list_observation_ids()[0]
        grids = vol_file.list_grids(obs_id)
        tensor_comps = vol_file.list_tensor_components(obs_id, grids[0])
        self.assertEqual(tensor_comps[4], 'Psi')
        tensor_data = vol_file.get_tensor_component(
            obs_id, grids[0], tensor_comps[4])
        target_tensor_data = DataStructures.DataVector([-0.0168989275087,
                                                        0.999906147623,
                                                        0.999906147623,
                                                        0.0171151930485,
                                                        0.999906147623,
                                                        0.0171151930485,
                                                        0.0171151930485,
                                                        -0.999905612058])
        self.assertEqual(tensor_data, target_tensor_data)
        h5_file.close()


if __name__ == '__main__':
    unittest.main(verbosity=2)
