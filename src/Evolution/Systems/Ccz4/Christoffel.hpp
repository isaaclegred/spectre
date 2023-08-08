// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Utilities/Gsl.hpp"

namespace Ccz4 {
/// @{
/*!
 * \brief Computes the conformal spatial christoffel symbols of the second kind.
 *
 * \details Computes the christoffel symbols as:
 * \f{align}
 *     \tilde{\Gamma}^k_{ij} &=
 *         \tilde{\gamma}^{kl} (D_{ijl} + D_{jil} - D_{lij})
 * \f}
 * where \f$\tilde{\gamma}^{ij}\f$ and \f$D_{ijk}\f$ are the inverse conformal
 * spatial metric and the CCZ4 auxiliary variable defined by
 * `Ccz4::Tags::InverseConformalMetric` and `Ccz4::Tags::FieldD`, respectively.
 */
template <typename DataType, size_t Dim, typename Frame>
void conformal_christoffel_second_kind(
    const gsl::not_null<tnsr::Ijj<DataType, Dim, Frame>*> result,
    const tnsr::II<DataType, Dim, Frame>& inverse_conformal_spatial_metric,
    const tnsr::ijj<DataType, Dim, Frame>& field_d);

template <typename DataType, size_t Dim, typename Frame>
tnsr::Ijj<DataType, Dim, Frame> conformal_christoffel_second_kind(
    const tnsr::II<DataType, Dim, Frame>& inverse_conformal_spatial_metric,
    const tnsr::ijj<DataType, Dim, Frame>& field_d);
/// @}

/// @{
/*!
 * \brief Computes the spatial christoffel symbols of the second kind.
 *
 * \details Computes the christoffel symbols as:
 * \f{align}
 *     \Gamma^k_{ij} &= \tilde{\Gamma}^k_{ij} -
 *         \tilde{\gamma}^{kl} (\tilde{\gamma}_{jl} P_i +
 *                              \tilde{\gamma}_{il} P_j -
 *                              \tilde{\gamma}_{ij} P_l)
 * \f}
 * where \f$\tilde{\gamma}^{ij}\f$, \f$\tilde{\gamma}_{ij}\f$,
 * \f$\tilde{\Gamma}^k_{ij}\f$, and \f$P_i\f$ are the conformal spatial metric,
 * the inverse conformal spatial metric, the conformal spatial christoffel
 * symbols of the second kind, and the CCZ4 auxiliary variable defined by
 * `Ccz4::Tags::ConformalMetric`, `Ccz4::Tags::InverseConformalMetric`,
 * `Ccz4::Tags::ConformalChristoffelSecondKind`, and `Ccz4::Tags::FieldP`,
 * respectively.
 */
template <typename DataType, size_t Dim, typename Frame>
void christoffel_second_kind(
    const gsl::not_null<tnsr::Ijj<DataType, Dim, Frame>*> result,
    const tnsr::ii<DataType, Dim, Frame>& conformal_spatial_metric,
    const tnsr::II<DataType, Dim, Frame>& inverse_conformal_spatial_metric,
    const tnsr::i<DataType, Dim, Frame>& field_p,
    const tnsr::Ijj<DataType, Dim, Frame>& conformal_christoffel_second_kind);

template <typename DataType, size_t Dim, typename Frame>
tnsr::Ijj<DataType, Dim, Frame> christoffel_second_kind(
    const tnsr::ii<DataType, Dim, Frame>& conformal_spatial_metric,
    const tnsr::II<DataType, Dim, Frame>& inverse_conformal_spatial_metric,
    const tnsr::i<DataType, Dim, Frame>& field_p,
    const tnsr::Ijj<DataType, Dim, Frame>& conformal_christoffel_second_kind);
/// @}

/// @{
/*!
 * \brief Computes the contraction of the conformal spatial Christoffel symbols
 * of the second kind.
 *
 * \details Computes the contraction as:
 *
 * \f{align}
 *     \tilde{\Gamma}^i &= \tilde{\gamma}^{jl} \tilde{\Gamma}^i_{jl}
 * \f}
 *
 * where \f$\tilde{\gamma}^{ij}\f$ is the inverse conformal spatial metric
 * defined by `Ccz4::Tags::InverseConformalMetric` and
 * \f$\tilde{\Gamma}^k_{ij}\f$ is the conformal spatial Christoffel symbols of
 * the second kind defined by `Ccz4::Tags::ConformalChristoffelSecondKind`.
 */
template <typename DataType, size_t Dim, typename Frame>
void contracted_conformal_christoffel_second_kind(
    const gsl::not_null<tnsr::I<DataType, Dim, Frame>*> result,
    const tnsr::II<DataType, Dim, Frame>& inverse_conformal_spatial_metric,
    const tnsr::Ijj<DataType, Dim, Frame>& conformal_christoffel_second_kind);

template <typename DataType, size_t Dim, typename Frame>
tnsr::I<DataType, Dim, Frame> contracted_conformal_christoffel_second_kind(
    const tnsr::II<DataType, Dim, Frame>& inverse_conformal_spatial_metric,
    const tnsr::Ijj<DataType, Dim, Frame>& conformal_christoffel_second_kind);
/// @}
}  // namespace Ccz4
