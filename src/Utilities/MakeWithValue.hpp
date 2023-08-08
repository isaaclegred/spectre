// Distributed under the MIT License.
// See LICENSE.txt for details.

///\file
/// Defines make_with_value

#pragma once

#include <array>
#include <complex>
#include <functional>
#include <type_traits>
#include <vector>

#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ForceInline.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \ingroup DataStructuresGroup
/// Implementations of make_with_value.
namespace MakeWithValueImpls {
/// Defines a method for determining the number of points represented
/// by an object.  This allows the object to appear as the input to
/// make_with_value.
///
/// The MakeWithValueImpls::number_of_points convenience wrapper is
/// provided to simplify calling this.
template <typename T, typename = std::nullptr_t>
struct NumberOfPoints {
  /// The default implementation will produce a compile-time error.
  [[noreturn]] static SPECTRE_ALWAYS_INLINE size_t apply(const T& /*input*/) {
    static_assert(typename tmpl::has_type<T, std::false_type>::type{},
                  "Do not know how to obtain a size from this type.  Either "
                  "implement NumberOfPoints or specialize MakeWithValueImpl "
                  "for the type you are trying to create.");
  }
};

/// The number of points represented by an object.
template <typename T>
size_t number_of_points(const T& input) {
  return NumberOfPoints<T>::apply(input);
}

/// Defines a method for producing an object representing a given
/// number of points.
///
/// Do not call these functions directly.  Use make_with_value
/// instead, which can take a size as its first argument.
template <typename R, typename = std::nullptr_t>
struct MakeWithSize {
  /// The default implementation will produce a compile-time error.
  /// In specializations, the \p value parameter need not be a template.
  template <typename T>
  [[noreturn]] static SPECTRE_ALWAYS_INLINE R apply(const size_t /*size*/,
                                                    const T& /*value*/) {
    static_assert(typename tmpl::has_type<R, std::false_type>::type{},
                  "Do not know how to create a sized object of this type.  "
                  "Either implement MakeWithSize or specialize "
                  "MakeWithValueImpl for the type you are trying to create.");
  }
};

template <typename R, typename T, typename = std::nullptr_t>
struct MakeWithValueImpl {
  /// The default implementation uses \ref number_of_points and MakeWithSize.
  template <typename ValueType>
  static SPECTRE_ALWAYS_INLINE R apply(const T& input, const ValueType value) {
    return MakeWithSize<R>::apply(number_of_points(input), value);
  }
};
}  // namespace MakeWithValueImpls

/// \ingroup DataStructuresGroup
/// \brief Given an object of type `T`, create an object of type `R` whose
/// elements are initialized to `value`.
///
/// \details This function is useful in function templates in order to
/// initialize the return type of a function template with `value` for functions
/// that can be called either at a single grid-point or to fill a data structure
/// at the same set of grid-points as the `input`

/// \tparam ValueType The type of `value`. For most containers, this will be
/// `double`.
///
/// \see MakeWithValueImpls, set_number_of_grid_points
template <typename R, typename T, typename ValueType>
SPECTRE_ALWAYS_INLINE std::remove_const_t<R> make_with_value(
    const T& input, const ValueType& value) {
  return MakeWithValueImpls::MakeWithValueImpl<std::remove_const_t<R>,
                                               T>::apply(input, value);
}

namespace MakeWithValueImpls {
template <>
struct NumberOfPoints<size_t> {
  static SPECTRE_ALWAYS_INLINE size_t apply(const size_t& input) {
    return input;
  }
};

/// \brief Returns a double initialized to `value` (`input` is ignored)
template <typename T>
struct MakeWithValueImpl<double, T> {
  static SPECTRE_ALWAYS_INLINE double apply(const T& /* input */,
                                            const double value) {
    return value;
  }
};

template <typename T>
struct MakeWithValueImpl<std::complex<double>, T> {
  static SPECTRE_ALWAYS_INLINE std::complex<double> apply(
      const T& /* input */, const std::complex<double> value) {
    return value;
  }
};

/// \brief Makes a `std::array`; each element of the `std::array`
/// must be `make_with_value`-creatable from a `InputType`.
template <size_t Size, typename T, typename InputType>
struct MakeWithValueImpl<std::array<T, Size>, InputType> {
  template <typename ValueType>
  static SPECTRE_ALWAYS_INLINE std::array<T, Size> apply(
      const InputType& input, const ValueType value) {
    return make_array<Size>(make_with_value<T>(input, value));
  }
};

template <size_t Size, typename T>
struct NumberOfPoints<std::array<T, Size>> {
  static SPECTRE_ALWAYS_INLINE size_t apply(const std::array<T, Size>& input) {
    static_assert(Size > 0);
    // size_t is interpreted as the number of points in other
    // contexts, but that doesn't make sense here.
    static_assert(not std::is_same_v<T, size_t>,
                  "Cannot get size from non-vector.");
    const size_t points = number_of_points(input[0]);
    ASSERT(
        alg::all_of(input,
                    [&](const T& t) { return number_of_points(t) == points; }),
        "Inconsistent number of points in array entries.");
    return points;
  }
};

template <typename T>
struct NumberOfPoints<std::vector<T>> {
  static SPECTRE_ALWAYS_INLINE size_t apply(const std::vector<T>& input) {
    // size_t is interpreted as the number of points in other
    // contexts, but that doesn't make sense here.
    static_assert(not std::is_same_v<T, size_t>,
                  "Cannot get number_of_points from non-vector.");
    ASSERT(not input.empty(),
           "Cannot get number of points from empty std::vector.");
    const size_t points = number_of_points(input[0]);
    ASSERT(
        alg::all_of(input,
                    [&](const T& t) { return number_of_points(t) == points; }),
        "Inconsistent number of points in vector entries.");
    return points;
  }
};

template <typename T>
struct NumberOfPoints<std::reference_wrapper<T>> {
  static SPECTRE_ALWAYS_INLINE size_t apply(
      const std::reference_wrapper<T>& input) {
    return number_of_points(input.get());
  }
};

/// \brief Makes a `TaggedTuple`; each element of the `TaggedTuple`
/// must be `make_with_value`-creatable from a `T`.
template <typename... Tags, typename T>
struct MakeWithValueImpl<tuples::TaggedTuple<Tags...>, T> {
  template <typename ValueType>
  static SPECTRE_ALWAYS_INLINE tuples::TaggedTuple<Tags...> apply(
      const T& input, const ValueType value) {
    return tuples::TaggedTuple<Tags...>(
        make_with_value<typename Tags::type>(input, value)...);
  }
};

template <typename Tag, typename... Tags>
struct NumberOfPoints<tuples::TaggedTuple<Tag, Tags...>> {
  static SPECTRE_ALWAYS_INLINE size_t apply(
      const tuples::TaggedTuple<Tag, Tags...>& input) {
    const size_t points = number_of_points(tuples::get<Tag>(input));
    ASSERT((... and (number_of_points(tuples::get<Tags>(input)) == points)),
           "Inconsistent number of points in tuple entries.");
    return points;
  }
};
}  // namespace MakeWithValueImpls
