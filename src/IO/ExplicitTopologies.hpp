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

Sphere<1>;
Sphere<2>;
Sphere<3>;

template <size_t Dim>
Torus public : TopologicalSpace{};

Torus<1>;
Torus<2>;
Torus<3>;

}  // namespace vis::detail
