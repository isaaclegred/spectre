// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Structure/Direction.hpp"

#include <ostream>
#include <pup.h>
#include <pup_stl.h>

#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"

template <size_t VolumeDim>
void Direction<VolumeDim>::pup(PUP::er& p) {
  size_t version = 0;
  p | version;
  // Remember to increment the version number when making changes to this
  // function. Retain support for unpacking data written by previous versions
  // whenever possible. See `Domain` docs for details.
  if (version >= 0) {
    p | axis_;
    p | side_;
  }
}

template <>
Direction<1>::Direction(const size_t dimension, const Side side) {
  ASSERT(
      0 == dimension,
      "dim = " << dimension << ", for Direction<1> only dim = 0 is allowed.");
  axis_ = Axis::Xi;
  side_ = side;
}

template <>
Direction<2>::Direction(const size_t dimension, const Side side) {
  ASSERT(0 == dimension or 1 == dimension,
         "dim = " << dimension
                  << ", for Direction<2> only dim = 0 or dim = 1 are allowed.");
  axis_ = 0 == dimension ? Axis::Xi : Axis::Eta;
  side_ = side;
}

template <>
Direction<3>::Direction(const size_t dimension, const Side side) {
  ASSERT(0 == dimension or 1 == dimension or 2 == dimension,
         "dim = " << dimension
                  << ", for Direction<3> only dim = 0, dim = 1, "
                     "or dim = 2 are allowed.");
  if (0 == dimension) {
    axis_ = Axis::Xi;
  }
  if (1 == dimension) {
    axis_ = Axis::Eta;
  }
  if (2 == dimension) {
    axis_ = Axis::Zeta;
  }
  side_ = side;
}

template <size_t VolumeDim>
std::ostream& operator<<(std::ostream& os,
                         const Direction<VolumeDim>& direction) {
  if (-1.0 == direction.sign()) {
    os << "-";
  } else {
    os << "+";
  }
  os << direction.dimension();
  return os;
}

template <size_t VolumeDim>
bool operator<(const Direction<VolumeDim>& lhs,
               const Direction<VolumeDim>& rhs) {
  if (lhs.axis() != rhs.axis()) {
    return lhs.axis() < rhs.axis();
  }
  return lhs.side() < rhs.side();
}

#define GET_DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data)                                        \
  template std::ostream& operator<<(std::ostream&,                    \
                                    const Direction<GET_DIM(data)>&); \
  template bool operator<(const Direction<GET_DIM(data)>& lhs,        \
                          const Direction<GET_DIM(data)>& rhs);       \
  template void Direction<GET_DIM(data)>::pup(PUP::er&);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef GET_DIM
#undef INSTANTIATION
