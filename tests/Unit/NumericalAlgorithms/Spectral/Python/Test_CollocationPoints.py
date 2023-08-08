# Distributed under the MIT License.
# See LICENSE.txt for details.

import unittest

import numpy as np
from numpy.polynomial import chebyshev, legendre

from spectre.DataStructures import DataVector
from spectre.Spectral import Basis, Mesh1D, Quadrature, collocation_points


class TestCollocationPoints(unittest.TestCase):
    def test_collocation_points_legendre_gauss(self):
        for num_coll_points in range(2, 13):
            coefs = np.zeros(num_coll_points + 1)
            coefs[-1] = 1.0
            coll_points_numpy = legendre.legroots(coefs)
            mesh = Mesh1D(num_coll_points, Basis.Legendre, Quadrature.Gauss)
            coll_points_spectre = collocation_points(mesh)
            np.testing.assert_allclose(
                coll_points_numpy, coll_points_spectre, 1e-12, 1e-12
            )

    def test_collocation_points_chebyshev_gauss(self):
        for num_coll_points in range(2, 13):
            coefs = np.zeros(num_coll_points + 1)
            coefs[-1] = 1.0
            coll_points_numpy = chebyshev.chebroots(coefs)
            mesh = Mesh1D(num_coll_points, Basis.Chebyshev, Quadrature.Gauss)
            coll_points_spectre = collocation_points(mesh)
            np.testing.assert_allclose(
                coll_points_numpy, coll_points_spectre, 1e-12, 1e-12
            )

    def test_collocation_points_legendre_gauss_lobatto(self):
        for num_coll_points in range(2, 13):
            coefs = np.zeros(num_coll_points)
            coefs[-1] = 1.0
            coll_points_numpy = np.concatenate(
                ([-1.0], legendre.Legendre(coefs).deriv().roots(), [1.0])
            )
            mesh = Mesh1D(
                num_coll_points, Basis.Legendre, Quadrature.GaussLobatto
            )
            coll_points_spectre = collocation_points(mesh)
            np.testing.assert_allclose(
                coll_points_numpy, coll_points_spectre, 1e-12, 1e-12
            )

    def test_collocation_points_chebyshev_gauss_lobatto(self):
        for num_coll_points in range(2, 13):
            coefs = np.zeros(num_coll_points)
            coefs[-1] = 1.0
            coll_points_numpy = np.concatenate(
                ([-1.0], chebyshev.Chebyshev(coefs).deriv().roots(), [1.0])
            )
            mesh = Mesh1D(
                num_coll_points, Basis.Chebyshev, Quadrature.GaussLobatto
            )
            coll_points_spectre = collocation_points(mesh)
            np.testing.assert_allclose(
                coll_points_numpy, coll_points_spectre, 1e-12, 1e-12
            )


if __name__ == "__main__":
    unittest.main()
