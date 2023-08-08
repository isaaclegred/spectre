// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"
#include "Utilities/TMPL.hpp"

/// @{
/*!
 * \ingroup TensorGroup
 * \brief Compute the Euclidean magnitude of a rank-1 tensor
 *
 * \details
 * Computes the square root of the sum of the squares of the components of
 * the rank-1 tensor.
 */
template <typename DataType, typename Index>
Scalar<DataType> magnitude(
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector) {
  return Scalar<DataType>{sqrt(get(dot_product(vector, vector)))};
}

template <typename DataType, typename Index>
void magnitude(const gsl::not_null<Scalar<DataType>*> magnitude,
               const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector) {
  set_number_of_grid_points(magnitude, vector);
  dot_product(magnitude, vector, vector);
  get(*magnitude) = sqrt(get(*magnitude));
}
/// @}

/// @{
/*!
 * \ingroup TensorGroup
 * \brief Compute the magnitude of a rank-1 tensor
 *
 * \details
 * Returns the square root of the input tensor contracted twice with the given
 * metric.
 */
template <typename DataType, typename Index>
Scalar<DataType> magnitude(
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector,
    const Tensor<DataType, Symmetry<1, 1>,
                 index_list<change_index_up_lo<Index>,
                            change_index_up_lo<Index>>>& metric) {
  Scalar<DataType> local_magnitude{get_size(get<0>(vector))};
  magnitude(make_not_null(&local_magnitude), vector, metric);
  return local_magnitude;
}

template <typename DataType, typename Index>
void magnitude(const gsl::not_null<Scalar<DataType>*> magnitude,
               const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector,
               const Tensor<DataType, Symmetry<1, 1>,
                            index_list<change_index_up_lo<Index>,
                                       change_index_up_lo<Index>>>& metric) {
  dot_product(magnitude, vector, vector, metric);
  get(*magnitude) = sqrt(get(*magnitude));
}
/// @}

namespace Tags {
/// \ingroup DataBoxTagsGroup
/// \ingroup DataStructuresGroup
/// The magnitude of a (co)vector
template <typename Tag>
struct Magnitude : db::PrefixTag, db::SimpleTag {
  using tag = Tag;
  using type = Scalar<DataVector>;
};

/// \ingroup DataBoxTagsGroup
/// \ingroup DataStructuresGroup
/// The Euclidean magnitude of a (co)vector
///
/// This tag inherits from `Tags::Magnitude<Tag>`
template <typename Tag>
struct EuclideanMagnitude : Magnitude<Tag>, db::ComputeTag {
  using base = Magnitude<Tag>;
  using return_type = typename base::type;
  static constexpr auto function =
      static_cast<void (*)(const gsl::not_null<return_type*>,
                           const typename Tag::type&)>(&magnitude);
  using argument_tags = tmpl::list<Tag>;
};

/// \ingroup DataBoxTagsGroup
/// \ingroup DataStructuresGroup
/// The magnitude of a (co)vector with respect to a specific metric
///
/// This tag inherits from `Tags::Magnitude<Tag>`
template <typename Tag, typename MetricTag>
struct NonEuclideanMagnitude : Magnitude<Tag>, db::ComputeTag {
  using base = Magnitude<Tag>;
  using return_type = typename base::type;
  static constexpr auto function = static_cast<void (*)(
      const gsl::not_null<return_type*>, const typename Tag::type&,
      const typename MetricTag::type&)>(&magnitude);
  using argument_tags = tmpl::list<Tag, MetricTag>;
};

/// \ingroup DataBoxTagsGroup
/// \ingroup DataStructuresGroup
/// The normalized (co)vector represented by Tag
template <typename Tag>
struct Normalized : db::PrefixTag, db::SimpleTag {
  using tag = Tag;
  using type = typename Tag::type;
};

/// \ingroup DataBoxTagsGroup
/// \ingroup DataStructuresGroup
/// Normalizes the (co)vector represented by Tag
///
/// This tag inherits from `Tags::Normalized<Tag>`
template <typename Tag>
struct NormalizedCompute : Normalized<Tag>, db::ComputeTag {
  using base = Normalized<Tag>;
  using return_type = typename base::type;
  static void function(const gsl::not_null<return_type*> normalized_vector,
                       const typename Tag::type& vector_in,
                       const typename Magnitude<Tag>::type& magnitude) {
    *normalized_vector = vector_in;
    for (size_t d = 0; d < normalized_vector->index_dim(0); ++d) {
      normalized_vector->get(d) /= get(magnitude);
    }
  }
  using argument_tags = tmpl::list<Tag, Magnitude<Tag>>;
};

}  // namespace Tags
