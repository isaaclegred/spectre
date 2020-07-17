// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>

#include "Connectivity.hpp"

namespace vis::detail {
template <size_t Dim>
Sphere public : TopologicalSpace{};
  
template<>
class Sphere<1>;

template<>
class Sphere<2>;

template<>
class Sphere<3>;

template <size_t Dim>
Torus public : TopologicalSpace{};

template<>
class Torus<1>;

template<>
class Torus<2>;

template<>
class Torus<3>;

}  // namespace vis::detail
