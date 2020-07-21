// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>

#include "Connectivity.hpp"

namespace vis::detail {

TopologicalSpace space_from_tag(Topology top, std::vector<size_t> extents);

template<size_t dim>
class Sphere : public  TopologicalSpace{
public:
  Sphere(std::vector<size_t> in_extents) : TopologicalSpace(in_extents){};
  std::vector<CellInTopology> topology_cells();

};

template 
class Sphere<1>;

template
class Sphere<2>;
  
template
class Sphere<3>;


template<size_t dim>
class Euclidean :  public TopologicalSpace{
public:
  Euclidean(std::vector<size_t> in_extents) : TopologicalSpace(in_extents){};
  std::vector<CellInTopology> topology_cells();
};

template
class Euclidean<1>;

template
class Euclidean<2>;

template
class Euclidean<3>;


}  // namespace vis::detail
