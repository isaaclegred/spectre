// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "IO/Connectivity.hpp"

#include <algorithm>

#include "DataStructures/Index.hpp"  // IWYU pragma: keep
#include "ErrorHandling/Error.hpp"

namespace vis {
namespace detail {
std::ostream& operator<<(std::ostream& os,
                         const BasicTopology& topology) noexcept {
  switch (topology) {
    case BasicTopology::Line:
      return os << "Line";
    case BasicTopology::Quad:
      return os << "Quad";
    case BasicTopology::Hexahedron:
      return os << "Hexahedron";
    // LCOV_EXCL_START
    default:
      ERROR("Unknown value of BasicTopology trying to be streamed");
      // LCOV_EXCL_STOP
  }
}

std::vector<CellInBasicTopology> cells_in_i1(
    const size_t number_of_points_in_i1) {
  // The number of cells in an I1 with N collocation points is N-1. Each cell
  // goes from the ith to the i+1th collocation point. This is stored in a
  // vector.
  std::vector<CellInBasicTopology> result;
  result.reserve(number_of_points_in_i1 - 1);
  for (size_t i = 0; i < number_of_points_in_i1 - 1; ++i) {
    result.emplace_back(BasicTopology::Line, std::vector<size_t>{i, i + 1});
  }
  return result;
}

std::vector<CellInBasicTopology> tensor_product_cells(
    const std::vector<CellInBasicTopology>& first_topology_cells,
    const size_t first_topology_cells_size,
    const std::vector<CellInBasicTopology>& second_topology_cells) {
  // first_topology_cells_size is the number of cells plus one in the topology.
  // For an I1 topology it is the number of grid points. Passing it in
  // separately allows easier generalization to other topologies such as S2.
  std::vector<CellInBasicTopology> result;
  for (const auto& first_cells : first_topology_cells) {
    for (const auto& second_cells : second_topology_cells) {
      const std::vector<size_t>& first_bounding_indices =
          first_cells.bounding_indices;
      const std::vector<size_t> second_bounding_indices =
          [first_topology_cells_size](
              const std::vector<size_t>& bounding_indices) {
            std::vector<size_t> bounding_computed(bounding_indices.size());
            std::transform(bounding_indices.begin(), bounding_indices.end(),
                           bounding_computed.begin(),
                           [first_topology_cells_size](const size_t element) {
                             return element * first_topology_cells_size;
                           });
            return bounding_computed;
          }(second_cells.bounding_indices);

      if (first_cells.topology == BasicTopology::Line and
          second_cells.topology == BasicTopology::Line) {
        result.emplace_back(
            BasicTopology::Quad,
            std::vector<size_t>{
                first_bounding_indices[0] + second_bounding_indices[0],
                first_bounding_indices[1] + second_bounding_indices[0],
                first_bounding_indices[1] + second_bounding_indices[1],
                first_bounding_indices[0] + second_bounding_indices[1]});
      } else if (first_cells.topology == BasicTopology::Line and
                 second_cells.topology == BasicTopology::Quad) {
        result.emplace_back(
            BasicTopology::Hexahedron,
            std::vector<size_t>{
                first_bounding_indices[0] + second_bounding_indices[0],
                first_bounding_indices[1] + second_bounding_indices[0],
                first_bounding_indices[1] + second_bounding_indices[1],
                first_bounding_indices[0] + second_bounding_indices[1],
                first_bounding_indices[0] + second_bounding_indices[3],
                first_bounding_indices[1] + second_bounding_indices[3],
                first_bounding_indices[1] + second_bounding_indices[2],
                first_bounding_indices[0] + second_bounding_indices[2]});
      } else {
        // LCOV_EXCL_START
        ERROR(
            "Tensor product of cells that are unknown. The first topology is: "
            << first_cells.topology
            << " and the second is: " << second_cells.topology);
        // LCOV_EXCL_STOP
      }
    }
  }
  return result;
}

template <size_t Dim>
std::vector<CellInBasicTopology> compute_cells(
    const Index<Dim>& extents) noexcept {
  std::vector<CellInBasicTopology> cells;
  std::vector<std::vector<CellInBasicTopology>> cells_per_topology;
  std::vector<size_t> size_per_topology;
  // Compute the number of cells in each topology and the number of extents
  // The extents are copied since this will need generalization for other
  // topologies.
  for (size_t i = 0; i < Dim; ++i) {
    cells_per_topology.push_back(cells_in_i1(extents[i]));
    size_per_topology.push_back(extents[i]);
  }
  cells = cells_per_topology.back();
  // Compute the tensor product of the number of cells in each topology to get
  // the cell structure of the full domain, ie 3D cell structure.
  for (size_t i = cells_per_topology.size() - 1; i > 0; --i) {
    cells = tensor_product_cells(cells_per_topology[i - 1],
                                 size_per_topology[i - 1], cells);
  }
  return cells;
}

std::vector<CellInBasicTopology> compute_cells(
    const std::vector<size_t>& extents) noexcept {
  if (extents.size() == 1) {
    return compute_cells(Index<1>{extents[0]});
  } else if (extents.size() == 2) {
    return compute_cells(Index<2>{extents[0], extents[1]});
  } else if (extents.size() == 3) {
    return compute_cells(Index<3>{extents[0], extents[1], extents[2]});
  }
  ERROR(
      "Only know how to compute connectivity for extents of size 1, 2, and 3, "
      "not "
      << extents.size());
}

// Explicit instantiations
template std::vector<CellInBasicTopology> compute_cells<1>(
    const Index<1>& extents) noexcept;
template std::vector<CellInBasicTopology> compute_cells<2>(
    const Index<2>& extents) noexcept;
template std::vector<CellInBasicTopology> compute_cells<3>(
    const Index<3>& extents) noexcept;

std::vector<CellInTopology> compute_cells(const TopologicalSpace& topology) noexcept{
  return topology.compute_topology();
  }


 

//Default implementation  
std::vector<CellInTopology> TopologicalSpace::compute_topology() const noexcept{
  return compute_cells(extents);
}
  
}  // namespace detail
}  // namespace vis
