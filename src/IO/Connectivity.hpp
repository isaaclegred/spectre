// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines functions for computing the connectivity of an element

#pragma once

#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>  // Distributed under the MIT License.
// See LICENSE.txt for details.

template <size_t Dim>
class Index;

/// Holds functions needed for visualizing data
namespace vis {
namespace detail {
/*!
 * \brief A list of all topologies for which we can compute the number of cells
 */
enum class BasicTopology { Line, Quad, Hexahedron };

enum class Topology {S1, S2, S3, E1, E2, E3};

std::ostream& operator<<(std::ostream& os, const BasicTopology& topology) noexcept;

/*!
 * \brief Represents the number of cells in a particular topology
 *
 * Each `CellInBasicTopology` holds an enum of type `BasicTopology` whose
 * value denotes the type of the topology, e.g. line, quad or hexahedron, and a
 * vector of bounding indices which are the indices of the grid coordinates in
 * the contiguous arrays of x, y, and z coordinates that bound the cell.
 */
struct CellInBasicTopology {
  // cppcheck-suppress passedByValue
  CellInBasicTopology(const BasicTopology& top,
                      std::vector<size_t> bounding_ind)
      : topology(top), bounding_indices(std::move(bounding_ind)) {}
  CellInBasicTopology() = default;
  CellInBasicTopology(const CellInBasicTopology& /*rhs*/) = default;
  CellInBasicTopology(CellInBasicTopology&& /*rhs*/) = default;
  CellInBasicTopology& operator=(const CellInBasicTopology& /*rhs*/) = default;
  CellInBasicTopology& operator=(CellInBasicTopology&& /*rhs*/) = default;
  ~CellInBasicTopology() = default;
  BasicTopology topology{BasicTopology::Line};
  std::vector<size_t> bounding_indices{};
};

using CellInTopology = CellInBasicTopology;
  
class TopologicalSpace {
public:
  TopologicalSpace(std::vector<size_t> in_extents) : extents(in_extents){};
  
  // Generalized notion of extents (I guess ultimately based on maps to R^n)
  std::vector<size_t> extents{};
  virtual ~TopologicalSpace() = default;
  virtual std::vector<CellInBasicTopology> compute_topology() const noexcept;
  
    
};

  
// @{
/*!
 * \brief Compute the cells in the element.
 *
 * Returns a vector of the cells in the topology I1^Dim, i.e. a line if Dim ==
 * 1, or a hexahedron if Dim == 3. The cells are bounded by lines connecting
 * grid points along the axes of the element, so if you have (n_x by n_y by n_z)
 * grid points, you have ((n_x-1) by (n_y-1) by (n_z-1)) cells.
 *
 * \note As more topologies are added, e.g. S2, the interface will need slight
 * modification, however the return type is likely to be able to remain the
 * same.
 */
template <size_t Dim>
std::vector<CellInBasicTopology> compute_cells(
    const Index<Dim>& extents) noexcept;

std::vector<CellInBasicTopology> compute_cells(
    const std::vector<size_t>& extents) noexcept;


std::vector<CellInTopology> compute_cells(const TopologicalSpace& topology) noexcept;

// @}
}  // namespace detail
}  // namespace vis
