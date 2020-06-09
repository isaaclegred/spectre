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
    def test_get_name(self):
        tc = TensorComponent("tensor component", DataVector(1.0, 1.1))
        self.assertEqual(tc.get_name(), "tensor component")

    def test_set_name(self):
        tc = TensorComponent("tensor component", DataVector(1.0, 1.1))
        tc.set_name("new tensor component")
        self.assertEqual(tc.get_name(), "new tensor component")

    def test_get_data(self):
        tc = TensorComponent("tensor component", DataVector(1.0, 1.1))
        self.assertEqual(tc.get_data(), [1.0, 1.1])

    def test_set_data(self):
        tc = TensorComponent("tensor component", DataVector(1.0, 1.1))
        tc.set_data(7.1, 5)
        self.assertEqual(tc.get_data(), [7.1, 5])

    def test_string(self):
        tc = TensorComponent("tensor component", DataVector(1.0, 1.1))
        self.assertEqual(str(tc), "(tensor component, (1.0, 1.1))")

    # Tests for ExtentsAndTensorVolumeData functions
    def test_get_extents(self):
        tc1 = TensorComponent("tensor component one", DataVector(1.0, 1.1))
        tc2 = TensorComponent("tensor component two", DataVector(7.1, 5))
        etvd = ExtentsAndTensorVolumeData([1, 2, 3, 4], [tc1, tc2])
        self.assertEqual(etvd.get_extents(), [1, 2, 3, 4])

    def test_set_extents(self):
        etvd = ExtentsAndTensorVolumeData([], [])
        etvd.set_extents([5, 6, 7, 8])
        self.assertEqual(etvd.get_extents(), [5, 6, 7, 8])

    def test_get_tensor_components(self):
        tc1 = TensorComponent("tensor component one", DataVector(1.0, 1.1))
        tc2 = TensorComponent("tensor component two", DataVector(7.1, 5))
        etvd = ExtentsAndTensorVolumeData([1, 2, 3, 4], [tc1, tc2])
        elf.assertEqual(etvd.get_tensor_components(), [tc1, tc2])

    def test_set_tensor_components(self):
        tc1 = TensorComponent("tensor component one", DataVector(1.0, 1.1))
        tc2 = TensorComponent("tensor component two", DataVector(7.1, 5))
        etvd = ExtentsAndTensorVolumeData([], [])
        etvd.set_tensor_components([tc1, tc2])
        self.assertEqual(etvd.get_tensor_components(), [tc1, tc2])
