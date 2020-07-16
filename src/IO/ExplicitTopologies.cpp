// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ExplicitTopologies.hpp"

#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>

#include "Connectivity.hpp"

namespace vis::detail {

std::vector<CellInTopolgy> Sphere<1>::compute_cells() {
  size_t num_points = extents[0];
  auto line_cells = compute_cells(num_points);
  // Attach the two loose ends of a line
  basic_cells.emplace_back(BasicTopology::Line,
                           std::vector<size_t>{num_points, 1});
  return basic_cells;
}

std::vector<CellInTopolgy> Sphere<2>::compute_cells() {
  size_t theta_pts = extents[0];
  size_t phi_pts = extents[1];
  auto square_cells = compute_cells(theta_pts, phi_pts);
  // Glue the theta points together with squares, and collapse the phi points
  // together at poles
}

}  // namespace vis::detail
