# Distributed under the MIT License.
# See LICENSE.txt for details.

from spectre.DataStructures import DataVector
from spectre.DataStructures import ExtentsAndTensorVolumeData
from spectre.DataStructures import TensorComponent
import unittest
import numpy as np
import numpy.testing as npt


class TestTensorData(unittest.TestCase):
    # Tests for TensorComponent functions
    def test_tensor_component(self):
        # Test name
        tc1 = TensorComponent("tensor component", DataVector([1.5, 1.1]))
        self.assertEqual(tc1.name, "tensor component")
        tc1.name = "new tensor component"
        self.assertEqual(tc1.name, "new tensor component")

        # Test data
        tc2 = TensorComponent("tensor component", DataVector([1.5, 1.1]))
        npt.assert_array_almost_equal(np.array(tc2.data), np.array([1.5, 1.1]))
        tc2.data = DataVector([7.1, 5])
        npt.assert_array_almost_equal(np.array(tc2.data), np.array([7.1, 5]))

        # Test string
        tc3 = TensorComponent("tensor component", DataVector([1.5, 1.1]))
        self.assertEqual(str(tc3), "(tensor component, (1.5,1.1))")

    # Tests for ExtentsAndTensorVolumeData functions
    def test_extents_and_tensor_components(self):
        # Test extents
        tc1 = TensorComponent("tensor component one", DataVector([1.5, 1.1]))
        tc2 = TensorComponent("tensor component two", DataVector([7.1, 5]))
        etvd1 = ExtentsAndTensorVolumeData([1, 2, 3, 4], [tc1, tc2])
        self.assertEqual(etvd1.extents, [1, 2, 3, 4])
        etvd1.extents = [5, 6, 7, 8]
        self.assertEqual(etvd1.extents, [5, 6, 7, 8])

        # Test tensor components
        tc3 = TensorComponent("tensor component three", DataVector([5.4, 8.2]))
        tc4 = TensorComponent("tensor component four", DataVector([0.9, 4.3]))
        etvd2 = ExtentsAndTensorVolumeData([1, 2, 3, 4], [tc1, tc2])
        expected_tensor_components1 = [tc1, tc2]
        for i in range(2):
            self.assertEqual(etvd2.tensor_components[i].name,
                             expected_tensor_components1[i].name)
            npt.assert_array_almost_equal(etvd2.tensor_components[i].data,
                             expected_tensor_components1[i].data)

        etvd2.tensor_components = [tc3, tc4]
        expected_tensor_components2 = [tc3, tc4]
        for i in range(2):
            self.assertEqual(etvd2.tensor_components[i].name,
                             expected_tensor_components2[i].name)
            npt.assert_array_almost_equal(etvd2.tensor_components[i].data,
                            expected_tensor_components2[i].data)


if __name__ == '__main__':
    unittest.main(verbosity=2)
