// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ExplicitTopologies.hpp"

#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>

#include "Connectivity.hpp"

namespace vis::detail {


template<>  
std::vector<CellInTopology> Sphere<1>::compute_topology() {
  size_t num_points = extents[0];
  if (num_points < 2) {
    ERROR("Constructing a 1-sphere requires at least 2 points");
  }
  auto line_cells = compute_cells(num_points);
  // Attach the two loose ends of a line
  basic_cells.emplace_back(BasicTopology::Line,
                           std::vector<size_t>{num_points, 1});
  return basic_cells;
}

template<>
std::vector<CellInTopology> Sphere<2>::compute_topology() {
  size_t theta_pts = extents[0];
  size_t phi_pts = extents[1];
  if (theta_pts < 2 or phi_pts < 2) {
    ERROR("2-sphere requires at least 2 theta points and 2 phi points");
  }
  std::vector<CellInTopology> result;
  //
  result.reserve((theta_pts - 1) * (phi_pts))
      // Two distinguished points
      size_t south_pole = 0;
  size_t north_pole = theta_pts * phi_pts - 1;  // is at least 3
  // Glue the theta points together with squares, and collapse the phi points
  // together at poles, first cells varies more slowly then second cells
  size_t points_so_far = 0;
  for (size_t theta_pt = 0; theta_pt < theta_pts; theta_pt++) {
    if (theta_pt = south_pole) {
      for (size_t phi_pt = 0; phi_pt < phi_pts; phi_pt++) {
        // Add polar wedges
        result.emplace_back(BasicTopology::Quad,
                            {south_pole, south_pole + phi_pt,
                             south_pole + phi_pt + 1, south_pole})
      }
      points_so_far += 1;
    }
    // The level below the north pole
    else if (theta_pt == theta_pts - 1) {
      // Add polar wedges
      result.emplace_back(BasicTopology::Quad,
                          {north_pole, north_pole - phi_pt,
                           north_pole - phi_pt + 1, north_pole})
    } else {  // Neither of the pole special cases
      for (size_t phi_pt = 0; phi_pt < phi_pts; theta_pt++) {
        // Add the square which is up and to the right of this point
        result.emplace_back(
            BasicTopology::Quad,
            {
                points_so_far + phi_pt,                // Bottom left
                points_so_far + phi_pt + 1,            // Bottom right
                points_so_far + phi_pt + phi_pts + 1,  // Upper right
                points_so_far + phi_pt + phi_pts        // Upper left
            })
      }
      points_so_far += phi_pts;
    }
  }
  return result;
}

template<>
std::vector<CellInTopology> Euclidean<1>::compute_topology(){
  return vis::detail::compute_cells(extents);
}

template<>
std::vector<CellInTopology> Euclidean<2>::compute_topology(){
  return vis::detail::compute_cells(extents);
}

template<>
std::vector<CellInTopology> Euclidean<3>::compute_topology(){
  return vis::detail::compute_cells(extents);}

TopologicalSpace space_from_tag(const Topology& top, const std::vector<size_t>&
                                extents) noexcept {
  TopologicalSpace* ptop;
  switch(top){
  case Topology::E1 : *ptop = Euclidean<1> (extents); break;
  case Topology::E2 : *ptop = Euclidean<2> (extents); break;
  case Topology::E3 : *ptop =  Euclidean<3> (extents); break;  
  case Topology::S1 : *ptop = Sphere<1> (extents); break;
  case Topology::S2 : *ptop = Sphere<2> (extents); break;
  case Topology::S3 : *ptop = Sphere<3> (extents); break;
  case default : ERROR("Topology not known");
  }
};


}  // namespace vis::detail
