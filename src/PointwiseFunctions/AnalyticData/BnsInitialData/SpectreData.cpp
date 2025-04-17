// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/BnsInitialData/SpectreData.hpp"

#include <Exporter.hpp>  // The SpEC Exporter
#include <memory>
#include <pup.h>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "IO/External/InterpolateFromSpec.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/Hydro/SpecificEnthalpy.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/System/ParallelInfo.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace BnsInitialData::AnalyticData {

template <typename DataType>
using background_tags = tmpl::list<
    hydro::Tags::RestMassDensity<DataVector>,
    gr::Tags::InverseSpatialMetric<DataType, 3>,
    gr::Tags::SpatialChristoffelSecondKindContracted<DataType, 3>,
    gr::Tags::Lapse<DataType>,
    ::Tags::deriv<gr::Tags::Lapse<DataType>, tmpl::integral_constant<size_t, 3>,
                  Frame::Inertial>,
    gr::Tags::Shift<DataType, 3>,
    ::Tags::deriv<gr::Tags::Shift<DataType, 3>,
                  tmpl::integral_constant<size_t, 3>, Frame::Inertial>,
    BnsInitialData::Tags::RotationalShift<DataType>,
    BnsInitialData::Tags::DerivLogLapseTimesDensityOverSpecificEnthalpy<
        DataType>,
    BnsInitialData::Tags::RotationalShiftStress<DataType>>;

template <size_t ThermodynamicDim>
SpectreData<ThermodynamicDim>::SpectreData(
    std::string volume_file_glob, std::string subfile_name,
    int observation_step,
    std::unique_ptr<equation_of_state_type> equation_of_state,
    const double density_cutoff, const double orbital_angular_velocity,
    const double euler_enthalpy_constant)
    : volume_file_glob_(std::move(volume_file_glob)),
      equation_of_state_(std::move(equation_of_state)),
      density_cutoff_(density_cutoff),
      orbital_angular_velocity_(orbital_angular_velocity),
      euler_enthalpy_constant_(euler_enthalpy_constant) {};

template <size_t ThermodynamicDim>
SpectreData<ThermodynamicDim>& SpectreData<ThermodynamicDim>::operator=(
    const SpectreData& rhs) {
  volume_file_glob_ = rhs.volume_file_glob_;
  subfile_name_ = rhs.subfile_name_;
  observation_step_ = rhs.observation_step_;
  equation_of_state_ = rhs.equation_of_state_->get_clone();
  density_cutoff_ = rhs.density_cutoff_;
  orbital_angular_velocity_ = rhs.orbital_angular_velocity_;
  euler_enthalpy_constant_ = rhs.euler_enthalpy_constant_;
  
  return *this;
}

template <size_t ThermodynamicDim>
SpectreData<ThermodynamicDim>::SpectreData(const SpectreData& rhs) {
  *this = rhs;
}

template <size_t ThermodynamicDim>
std::unique_ptr<elliptic::analytic_data::Background>
SpectreData<ThermodynamicDim>::get_clone() const {
  return std::make_unique<SpectreData>(*this);
}

template <size_t ThermodynamicDim>
SpectreData<ThermodynamicDim>::SpectreData(CkMigrateMessage* msg)
    : elliptic::analytic_data::Background(msg) {}

template <size_t ThermodynamicDim>
void SpectreData<ThermodynamicDim>::pup(PUP::er& p) {
  p | volume_file_glob_;
  p | subfile_name_;
  p | observation_step_;
  p | equation_of_state_;
  p | density_cutoff_;
  p | orbital_angular_velocity_;
  p | euler_enthalpy_constant_;
}

template <size_t ThermodynamicDim>
PUP::able::PUP_ID SpectreData<ThermodynamicDim>::my_PUP_ID = 0;

template <size_t ThermodynamicDim>
template <typename DataType>
tuples::tagged_tuple_from_typelist<
    typename SpectreData<ThermodynamicDim>::template interpolated_tags<DataType>>
SpectreData<ThermodynamicDim>::interpolate_from_spectre(
    const tnsr::I<DataType, 3>& x) const {
        return spectre::Exporter::interpolate_to_points<interpolated_tags>(
            volume_file_glob_,  subfile_name_, observation_step_, x, false, static_cast<size_t>(sys::my_local_rank()));
      }

// Deriv of velocity potential is only
// used for validation
template <size_t ThermodynamicDim>
template <typename DataType>
tnsr::i<DataType, 3> SpectreData<ThermodynamicDim>::deriv_of_velocity_potential(
    const tnsr::I<DataType, 3, Frame::Inertial>& x) const {
  const auto interpolated_vars = interpolate_from_spectre(x);
  const auto& lower_spatial_four_velocity =
      get<hydro::Tags::LowerSpatialFourVelocity<DataType, 3>>(
          interpolated_vars);
  auto rest_mass_density =
      get<hydro::Tags::RestMassDensity<DataType>>(interpolated_vars);
  get(rest_mass_density) += 1e-16;
  const DataType specific_enthalpy = select(
      step_function(get(rest_mass_density) - 0),
      get(equation_of_state_->pressure_from_density(rest_mass_density)) /
              get(rest_mass_density) +
          (1.0 + get(equation_of_state_->specific_internal_energy_from_density(
                     rest_mass_density))),
      make_with_value<DataType>(
          get(rest_mass_density),
          equation_of_state_->specific_enthalpy_lower_bound()));
  auto result = make_with_value<tnsr::i<DataType, 3>>(x, 0.0);

  tenex::evaluate<ti::i>(make_not_null(&result),
                         Scalar<DataVector>{specific_enthalpy}() *
                             lower_spatial_four_velocity(ti::i));
  return result;
}

// The velocity potential is used in the initial guess
template <size_t ThermodynamicDim>
template <typename DataType>
tuples::TaggedTuple<Tags::VelocityPotential<DataType>>
SpectreData<ThermodynamicDim>::variables(
    const tnsr::I<DataType, 3, Frame::Inertial>& x,
    tmpl::list<Tags::VelocityPotential<DataType>> /*meta*/) const {
  // return velocity potential (only a guess)
  Scalar<DataType> velocity_potential =
      make_with_value<Scalar<DataType>>(x, orbital_angular_velocity_);
  // This is not a good guess, but it's better than nothing.
  // It's not clear to me this should actually be used
  get(velocity_potential) *= (get<1>(x) * get<0>(x));

  return {std::move(velocity_potential)};
}
// The fixed sources are used in initialization
template <size_t ThermodynamicDim>
template <typename DataType>
tuples::TaggedTuple<::Tags::FixedSource<Tags::VelocityPotential<DataType>>>
SpectreData<ThermodynamicDim>::variables(
    const tnsr::I<DataType, 3, Frame::Inertial>& x, const Mesh<3>& mesh,
    const InverseJacobian<DataType, 3, Frame::ElementLogical, Frame::Inertial>&
        inv_jacobian,
    tmpl::list<::Tags::FixedSource<Tags::VelocityPotential<DataType>>> /*meta*/)
    const {
  auto background_values =
      variables(x, mesh, inv_jacobian, background_tags<DataType>{});
  tuples::TaggedTuple<::Tags::FixedSource<Tags::VelocityPotential<DataType>>>
      result{};
  const auto& lapse = get<gr::Tags::Lapse<DataType>>(background_values);
  const auto& shift = get<gr::Tags::Shift<DataType, 3>>(background_values);
  const auto& rotational_shift =
      get<Tags::RotationalShift<DataType>>(background_values);
  const auto& deriv_log_lapse_times_density_over_specific_enthalpy =
      get<Tags::DerivLogLapseTimesDensityOverSpecificEnthalpy<DataType>>(
          background_values);

  const auto& deriv_of_lapse =
      get<::Tags::deriv<gr::Tags::Lapse<DataType>,
                        tmpl::integral_constant<size_t, 3>, Frame::Inertial>>(
          background_values);

  const auto& deriv_of_shift =
      get<::Tags::deriv<gr::Tags::Shift<DataType, 3>,
                        tmpl::integral_constant<size_t, 3>, Frame::Inertial>>(
          background_values);
  const auto& spatial_christoffel_second_kind_contracted =
      get<gr::Tags::SpatialChristoffelSecondKindContracted<DataType, 3>>(
          background_values);

  ::tenex::evaluate<>(
      make_not_null(
          &get<::Tags::FixedSource<Tags::VelocityPotential<DataType>>>(result)),
      -euler_enthalpy_constant_ *
          (1.0 / square(lapse()) * rotational_shift(ti::I) *
               deriv_log_lapse_times_density_over_specific_enthalpy(ti::i) -
           2.0 / cube(lapse()) * rotational_shift(ti::I) *
               deriv_of_lapse(ti::i) +
           1.0 / square(lapse()) * deriv_of_shift(ti::i, ti::I) +
           // Christoffel terms, assume the spatial rotational
           // killing vector is (spatially) covariantly constant
           1.0 / square(lapse()) *
               (shift(ti::I) *
                spatial_christoffel_second_kind_contracted(ti::i))));

  return result;
}
template <size_t ThermodynamicDim>
template <typename DataType>
tuples::TaggedTuple<gr::Tags::InverseSpatialMetric<DataType, 3>>
SpectreData<ThermodynamicDim>::variables(
    const tnsr::I<DataType, 3, Frame::Inertial>& x,
    tmpl::list<gr::Tags::InverseSpatialMetric<DataType, 3>> /*meta*/) const {
  // interpolate from spec, then set gamma
  const auto interpolated_vars = interpolate_from_spectre(x);

  const auto& spatial_metric =
      get<gr::Tags::SpatialMetric<DataType, 3>>(interpolated_vars);
  tuples::TaggedTuple<gr::Tags::InverseSpatialMetric<DataType, 3>> result{};
  get<gr::Tags::InverseSpatialMetric<DataType, 3>>(result) =
      determinant_and_inverse(spatial_metric).second;
  return result;
}
template <size_t ThermodynamicDim>
template <typename DataType>
tuples::tagged_tuple_from_typelist<background_tags<DataType>>
SpectreData<ThermodynamicDim>::variables(
    const tnsr::I<DataType, 3, Frame::Inertial>& x, const Mesh<3>& mesh,
    const InverseJacobian<DataType, 3, Frame::ElementLogical, Frame::Inertial>&
        inv_jacobian,
    background_tags<DataType> /*meta*/) const {
  // interpolate from spec, take num derivatives, return
  // Shift, lapse spatial metric imported
  auto result = tuples::tagged_tuple_from_typelist<background_tags<DataType>>{};
  const auto interpolated_vars = interpolate_from_spectre(x);
  const auto& spatial_metric =
      get<gr::Tags::SpatialMetric<DataType, 3>>(interpolated_vars);
  const auto spatial_metric_determinant_and_inverse =
      determinant_and_inverse(spatial_metric);
  const auto& inv_spatial_metric =
      spatial_metric_determinant_and_inverse.second;
  get<gr::Tags::InverseSpatialMetric<DataType, 3>>(result) = inv_spatial_metric;
  const auto sqrt_det_spatial_metric =
      Scalar<DataType>{sqrt(get(spatial_metric_determinant_and_inverse.first))};
  const auto deriv_sqrt_det_spatial_metric =
      partial_derivative(sqrt_det_spatial_metric, mesh, inv_jacobian);
  // Get the one contracted Christoffel needed for fluxes
  const auto spatial_christoffel_second_kind_contracted =
      tenex::evaluate<ti::i>(deriv_sqrt_det_spatial_metric(ti::i) /
                             sqrt_det_spatial_metric());
  get<gr::Tags::SpatialChristoffelSecondKindContracted<DataType, 3>>(result) =
      spatial_christoffel_second_kind_contracted;
  get<gr::Tags::Lapse<DataType>>(result) =
      get<gr::Tags::Lapse<DataType>>(interpolated_vars);

  // Get Lapse and shift derivatives
  get<::Tags::deriv<gr::Tags::Lapse<DataType>,
                    tmpl::integral_constant<size_t, 3>, Frame::Inertial>>(
      result) =
      partial_derivative(get<gr::Tags::Lapse<DataType>>(interpolated_vars),
                         mesh, inv_jacobian);
  get<gr::Tags::Shift<DataType, 3>>(result) =
      get<gr::Tags::Shift<DataType, 3>>(interpolated_vars);
  get<::Tags::deriv<gr::Tags::Shift<DataType, 3>,
                    tmpl::integral_constant<size_t, 3>, Frame::Inertial>>(
      result) =
      partial_derivative(get<gr::Tags::Shift<DataType, 3>>(interpolated_vars),
                         mesh, inv_jacobian);
  // Get the rotational shift + deriv of log lapse over enthalpy + stress
  const auto spatial_rotational_killing_vector =
      hydro::initial_data::irrotational_bns::spatial_rotational_killing_vector(
          x, orbital_angular_velocity_);
  const auto rotational_shift =
      hydro::initial_data::irrotational_bns::rotational_shift(
          get<gr::Tags::Shift<DataType, 3>>(interpolated_vars),
          spatial_rotational_killing_vector);
  get<Tags::RotationalShift<DataType>>(result) = rotational_shift;
  auto rest_mass_density =
      get<hydro::Tags::RestMassDensity<DataType>>(interpolated_vars);
  get(rest_mass_density) += 1.0e-16;
  get<hydro::Tags::RestMassDensity<DataType>>(result) = rest_mass_density;
  // The SpEC solution should have e.g. B&S Eq. 15.76 satisfied; however,
  // we do not assume it is satisfied.  We take the SpEC density
  // (equivalently enthalpy) profile to be fixed and compute the
  // velocity potential from the matter and spacetime profiles.
  const DataType enthalpy_density = select(
      step_function(get(rest_mass_density) - 0.0),
      get(equation_of_state_->pressure_from_density(rest_mass_density)) +
          get(rest_mass_density) *
              (1.0 +
               get(equation_of_state_->specific_internal_energy_from_density(
                   rest_mass_density))),
      make_with_value<DataType>(
          get(rest_mass_density),
          equation_of_state_->specific_enthalpy_lower_bound()));
  const auto deriv_log_lapse_times_density_over_specific_enthalpy =
      partial_derivative(
          Scalar<DataType>{
              log(get(get<gr::Tags::Lapse<DataType>>(interpolated_vars)) *
                  square(get(rest_mass_density)) / enthalpy_density)},
          mesh, inv_jacobian);
  get<Tags::DerivLogLapseTimesDensityOverSpecificEnthalpy<DataType>>(result) =
      deriv_log_lapse_times_density_over_specific_enthalpy;
  const auto rotational_shift_stress =
      hydro::initial_data::irrotational_bns::rotational_shift_stress(
          rotational_shift, get<gr::Tags::Lapse<DataType>>(interpolated_vars));
  get<Tags::RotationalShiftStress<DataType>>(result) = rotational_shift_stress;
  return result;
}

#define THERMODIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data)                                                \
  template class SpectreData<THERMODIM(data)>;                                   \
  template tuples::tagged_tuple_from_typelist<                                \
      typename SpectreData<THERMODIM(data)>::template interpolated_tags<double>> \
  SpectreData<THERMODIM(data)>::interpolate_from_spectre(                           \
      const tnsr::I<double, 3>& x) const;                                     \
  template tuples::tagged_tuple_from_typelist<typename SpectreData<THERMODIM(    \
      data)>::template interpolated_tags<DataVector>>                         \
  SpectreData<THERMODIM(data)>::interpolate_from_spectre(                           \
      const tnsr::I<DataVector, 3>& x) const;                                 \
  template tnsr::i<DataVector, 3>                                             \
  SpectreData<THERMODIM(data)>::deriv_of_velocity_potential(                     \
      const tnsr::I<DataVector, 3, Frame::Inertial>& x) const;                \
  template tuples::TaggedTuple<                                               \
      ::Tags::FixedSource<Tags::VelocityPotential<DataVector>>>               \
  SpectreData<THERMODIM(data)>::variables(                                       \
      const tnsr::I<DataVector, 3, Frame::Inertial>& x, const Mesh<3>& mesh,  \
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,             \
                            Frame::Inertial>& inv_jacobian,                   \
      tmpl::list<                                                             \
          ::Tags::FixedSource<Tags::VelocityPotential<DataVector>>> /*meta*/) \
      const;                                                                  \
  template tuples::TaggedTuple<gr::Tags::InverseSpatialMetric<DataVector, 3>> \
  SpectreData<THERMODIM(data)>::variables(                                       \
      const tnsr::I<DataVector, 3, Frame::Inertial>& x,                       \
      tmpl::list<gr::Tags::InverseSpatialMetric<DataVector, 3>> /*meta*/)     \
      const;                                                                  \
  template tuples::TaggedTuple<Tags::VelocityPotential<DataVector>>           \
  SpectreData<THERMODIM(data)>::variables(                                       \
      const tnsr::I<DataVector, 3, Frame::Inertial>& x,                       \
      tmpl::list<Tags::VelocityPotential<DataVector>> /*meta*/) const;        \
  template tuples::tagged_tuple_from_typelist<                                \
      SpectreData<THERMODIM(data)>::background_tags<DataVector>>                 \
  SpectreData<THERMODIM(data)>::variables(                                       \
      const tnsr::I<DataVector, 3, Frame::Inertial>& x, const Mesh<3>& mesh,  \
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,             \
                            Frame::Inertial>& inv_jacobian,                   \
      background_tags<DataVector> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1))

#undef INSTANTIATION
#undef THERMODIM
}  // namespace BnsInitialData::AnalyticData
