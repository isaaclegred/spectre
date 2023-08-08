// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \ingroup GeneralRelativityGroup
/// Holds functions related to general relativity.
namespace gr {
/// @{
/*!
 * \ingroup GeneralRelativityGroup
 * \brief Compute lapse from shift and spacetime metric
 *
 * \details Computes
 * \f{align}
 *    \alpha &= \sqrt{\beta^i g_{it}-g_{tt}}
 * \f}
 * where \f$ \alpha \f$, \f$ \beta^i\f$, and \f$g_{ab}\f$ are the lapse, shift,
 * and spacetime metric.
 * This can be derived, e.g., from Eqs. 2.121--2.122 of Baumgarte & Shapiro.
 */
template <typename DataType, size_t SpatialDim, typename Frame>
Scalar<DataType> lapse(
    const tnsr::I<DataType, SpatialDim, Frame>& shift,
    const tnsr::aa<DataType, SpatialDim, Frame>& spacetime_metric);

template <typename DataType, size_t SpatialDim, typename Frame>
void lapse(gsl::not_null<Scalar<DataType>*> lapse,
           const tnsr::I<DataType, SpatialDim, Frame>& shift,
           const tnsr::aa<DataType, SpatialDim, Frame>& spacetime_metric);
/// @}

namespace Tags {
/*!
 * \brief Compute item for lapse \f$\alpha\f$ from the spacetime metric
 * \f$g_{ab}\f$ and the shift \f$\beta^i\f$.
 *
 * \details Can be retrieved using `gr::Tags::Lapse`.
 */
template <typename DataType, size_t SpatialDim, typename Frame>
struct LapseCompute : Lapse<DataType>, db::ComputeTag {
  using argument_tags =
      tmpl::list<Shift<DataType, SpatialDim, Frame>,
                 SpacetimeMetric<DataType, SpatialDim, Frame>>;

  using return_type = Scalar<DataType>;

  static constexpr auto function =
      static_cast<void (*)(const gsl::not_null<Scalar<DataType>*> lapse,
                           const tnsr::I<DataType, SpatialDim, Frame>&,
                           const tnsr::aa<DataType, SpatialDim, Frame>&)>(
          &lapse<DataType, SpatialDim, Frame>);

  using base = Lapse<DataType>;
};
}  // namespace Tags
}  // namespace gr
