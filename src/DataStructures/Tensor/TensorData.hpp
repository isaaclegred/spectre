// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataVector.hpp"
// Could I forward Declare instead?
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

/*!
 * \ingroup DataStructuresGroup
 * \brief An untyped tensor component with a name for observation.
 *
 * The name should be a path inside an H5 file, typically starting with the name
 * of the volume subfile. For example,
 * `element_volume_data.vol/ObservationId[ID]/[ElementIdName]/psi_xx`.
 */
struct TensorComponent {
  TensorComponent() = default;
  TensorComponent(std::string n, DataVector d) noexcept
      : name(std::move(n)), data(std::move(d)) {}
  void pup(PUP::er& p) noexcept;  // NOLINT
  std::string name{};
  DataVector data{};
};

std::ostream& operator<<(std::ostream& os, const TensorComponent& t) noexcept;

bool operator==(const TensorComponent& lhs,
                const TensorComponent& rhs) noexcept;

bool operator!=(const TensorComponent& lhs,
                const TensorComponent& rhs) noexcept;

/*!
 * \ingroup DataStructuresGroup
 * \brief Holds the extents of the mesh and the tensor components on the mesh.
 *
 * The extents is a `std::vector<size_t>` where each element is the number of
 * grid points in the given dimension. The `TensorComponent`s must live on the
 * grid of the size of the extents. We use runtime extents instead of the
 * `Index` class because observers may write 1D, 2D, or 3D data in a 3D
 * simulation.
 */
struct ExtentsAndTensorVolumeData {
  ExtentsAndTensorVolumeData() = default;
  ExtentsAndTensorVolumeData(std::vector<size_t> exts,
                             std::vector<TensorComponent> components) noexcept
      : extents(std::move(exts)), tensor_components(std::move(components)) {}
  void pup(PUP::er& p) noexcept;  // NOLINT
  std::vector<size_t> extents{};
  std::vector<TensorComponent> tensor_components{};
};

// We need a datastructure to store the quadrature and the spectral basis as it
// is passed to the volume writer.
struct ElementVolumeData : ExtentsAndTensorVolumeData {
  ElementVolumeData() = default;
  ElementVolumeData(std::vector<size_t> extents,
                    std::vector<TensorComponent> components,
                    std::vector <Spectral::Basis> basis,
                    std::vector<Spectral::Quadrature> quadrature) noexcept
  :ExtentsAndTensorVolumeData(extents, components),
    basis(std::move(basis)), quadrature(std::move(quadrature)) {}
  void pup(PUP::er& p) noexcept;
  std::vector<size_t> extents{};
  std::vector<TensorComponent> tensor_components{};
  std::vector<Spectral::Basis> basis{};
  std::vector<Spectral::Quadrature> quadrature{};

};
