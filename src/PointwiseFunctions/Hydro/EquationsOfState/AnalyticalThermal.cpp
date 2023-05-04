// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/Hydro/EquationsOfState/AnalyticalThermal.hpp"

#include <cmath>
#include <utility>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/RootFinding/TOMS748.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/Factory.hpp"
#include "Utilities/ConstantExpressions.hpp"
namespace {
// Find positive root of a polynomial, assuming there is exactly one
double root_polynomial_positive(const std::vector<double>& coefficients,
                                double lower_bound, double upper_bound) {
  auto miss = [&coefficients](double x) {
    return evaluate_polynomial(coefficients, x);
  };
  const auto root_from_lambda =
      RootFinder::toms748(miss, lower_bound, upper_bound, 1.0e-14, 1.0e-15);
  return root_from_lambda;
}
DataVector root_polynomial_positive(const std::vector<DataVector>& coefficients,
                                    double lower_bound, double upper_bound) {
  DataVector solutions = make_with_value<DataVector>(coefficients[0], 0.0);
  for (size_t i = 0; i < solutions.size(); i++) {
    std::vector<double> pointwise_coefficients(coefficients.size());
    std::transform(coefficients.begin(), coefficients.end(),
                   pointwise_coefficients.begin(),
                   [&i](const auto& vec) { return vec[i]; });
    solutions[i] = root_polynomial_positive(pointwise_coefficients, lower_bound,
                                            upper_bound);
  }
  return solutions;
}
}  // namespace
namespace EquationsOfState {
template <typename ColdEquationOfState>
AnalyticalThermal<ColdEquationOfState>::AnalyticalThermal(
    ColdEquationOfState cold_eos, const double S0, const double L,
    const double gamma, const double n0, const double alpha)
    : cold_eos_(std::move(cold_eos)),
      S0_(S0),
      L_(L),
      gamma_(gamma),
      n0_(n0),
      alpha_(alpha) {
  eta_ = get_eta();
}

EQUATION_OF_STATE_MEMBER_DEFINITIONS(template <typename ColdEquationOfState>,
                                     AnalyticalThermal<ColdEquationOfState>,
                                     double, 3)
EQUATION_OF_STATE_MEMBER_DEFINITIONS(template <typename ColdEquationOfState>,
                                     AnalyticalThermal<ColdEquationOfState>,
                                     DataVector, 3)

template <typename ColdEquationOfState>
std::unique_ptr<
    EquationOfState<AnalyticalThermal<ColdEquationOfState>::is_relativistic, 3>>
AnalyticalThermal<ColdEquationOfState>::get_clone() const {
  auto clone = std::make_unique<AnalyticalThermal<ColdEquationOfState>>(*this);
  return std::unique_ptr<EquationOfState<is_relativistic, 3>>(std::move(clone));
}

template <typename ColdEquationOfState>
bool AnalyticalThermal<ColdEquationOfState>::operator==(
    const AnalyticalThermal<ColdEquationOfState>& rhs) const {
  return cold_eos_ == rhs.cold_eos_ and S0_ == rhs.S0_ and L_ == rhs.L_ and
         gamma_ == rhs.gamma_ and n0_ == rhs.n0_ and alpha_ == rhs.alpha_;
}
template <typename ColdEquationOfState>
bool AnalyticalThermal<ColdEquationOfState>::operator!=(
    const AnalyticalThermal<ColdEquationOfState>& rhs) const {
  return not(*this == rhs);
}

template <typename ColdEquationOfState>
bool AnalyticalThermal<ColdEquationOfState>::is_equal(
    const EquationOfState<is_relativistic, 3>& rhs) const {
  const auto& derived_ptr =
      dynamic_cast<const AnalyticalThermal<ColdEquationOfState>* const>(&rhs);
  return derived_ptr != nullptr and *derived_ptr == *this;
}

template <typename ColdEquationOfState>
AnalyticalThermal<ColdEquationOfState>::AnalyticalThermal(CkMigrateMessage* msg)
    : EquationOfState<is_relativistic, 3>(msg) {}

template <typename ColdEquationOfState>
void AnalyticalThermal<ColdEquationOfState>::pup(PUP::er& p) {
  EquationOfState<is_relativistic, 3>::pup(p);
  p | cold_eos_;
  p | S0_;
  p | L_;
  p | gamma_;
  p | n0_;
  p | alpha_;
  p | eta_;
}
template <class ColdEos>
double AnalyticalThermal<ColdEos>::get_eta() const {
  return 5.0 / 9.0 * (L_ - 3 * S0_ * gamma_) /
         ((cbrt(0.25) - 1) * (2.0 / 3.0 - gamma_) *
          baryonic_fermi_internal_energy(saturation_density_));
}
template <class ColdEos>
template <class DataType>
DataType AnalyticalThermal<ColdEos>::baryonic_fermi_internal_energy(
    const DataType& rest_mass_density) const {
  return square(hbar_over_baryon_mass_to_four_thirds_) / 2.0 *
         pow(3 * square(M_PI) * rest_mass_density, 2.0 / 3.0);
}

template <class ColdEos>
template <class DataType>
DataType AnalyticalThermal<ColdEos>::thermal_internal_energy(
    const DataType& rest_mass_density, const DataType& temperature,
    const DataType& electron_fraction) const {
  // Following box 1
  // TODO: Do the right thing here pointwise
  const double fs = max(temperature) < .001 ? 1.0 : 11.0 / 4.0;
  const DataType radiation = 4.0 * stefan_boltzmann_sigma_ * fs *
                             pow(temperature, 4) / rest_mass_density;
  // Factor out a `temperature` from the ideal and degenerate terms
  const double ideal = 1.5;

  const DataType degenerate =
      (a_degeneracy(rest_mass_density,
                    make_with_value<DataType>(rest_mass_density, 0.5),
                    dirac_effective_mass(rest_mass_density)) +
       a_degeneracy(rest_mass_density, electron_fraction,
                    make_with_value<DataType>(
                        rest_mass_density, electron_mass_over_baryon_mass_))) *
      temperature;
  return radiation + temperature * ideal * degenerate / (ideal + degenerate);
}
template <class ColdEos>
template <class DataType>
DataType AnalyticalThermal<ColdEos>::thermal_pressure(
    const DataType& rest_mass_density, const DataType& temperature,
    const DataType& electron_fraction) const {
  // Following box 2
  const double fs = max(temperature) < .001 ? 1.0 : 11.0 / 4.0;
  const DataType radiation =
      4.0 / 3.0 * stefan_boltzmann_sigma_ * fs * pow(temperature, 4);
  const DataType ideal = rest_mass_density;
  const DataType degenerate = rest_mass_density * temperature *
                              (a_degeneracy_log_density_derivative(
                                  rest_mass_density, electron_fraction,
                                  dirac_effective_mass(rest_mass_density) +
                                      a_degeneracy_log_density_derivative(
                                          rest_mass_density, electron_fraction,
                                          electron_mass_over_baryon_mass_)));
  return radiation + temperature * ideal * degenerate / (ideal + degenerate);
}
template <class ColdEos>
template <class DataType>
DataType AnalyticalThermal<ColdEos>::thermal_pressure_density_derivative(
    const DataType& rest_mass_density, const DataType& temperature,
    const DataType& electron_fraction) const {
  // Following Eq. (61)
  const double fs = max(temperature) < .001 ? 1.0 : 11.0 / 4.0;
  const DataType radiation = 16.0 * stefan_boltzmann_sigma_ * fs *
                             pow(temperature, 4) / (9.0 * rest_mass_density);
  // Factor out a T from the ideal and degenerate contributions
  const double ideal = 1.6;
  const DataType baryon_effective_mass =
      dirac_effective_mass(rest_mass_density);
  const DataType a_electron = a_degeneracy(rest_mass_density, electron_fraction,
                                           electron_mass_over_baryon_mass_);
  const DataType a_SM =
      a_degeneracy(rest_mass_density, electron_fraction, baryon_effective_mass);
  const DataType a_electron_log_density_deriv =
      a_degeneracy_log_density_derivative(rest_mass_density, electron_fraction,
                                          electron_mass_over_baryon_mass_);
  const DataType a_SM_log_density_deriv = a_degeneracy_log_density_derivative(
      rest_mass_density, electron_fraction, baryon_effective_mass);
  const DataType degenerate_temperature_dependent =
      2 * thermal_pressure(rest_mass_density, temperature, electron_fraction) /
      rest_mass_density *
      (1.0 * (a_SM_log_density_deriv + a_electron_log_density_deriv) *
       electron_fraction) /
      (a_SM + a_electron * electron_fraction);
  const DataType A = -alpha_ * (1 - square(baryon_effective_mass));

  const DataType C =
      pow(3.0 * square(M_PI) * electron_fraction * rest_mass_density,
          2.0 / 3.0) *
      square(hbar_over_baryon_mass_to_four_thirds_ / baryon_effective_mass);
  const DataType B = 1 / (1 + C);
  const DataType rho_times_a_SM_second_log_density_deriv =
      a_SM_log_density_deriv * (a_SM_log_density_deriv / a_SM) +
      2.0 * a_SM / 3.0 * B *
          (3.0 * square(A) - 1.0 / 3.0 * B * square(3.0 * A + C) +
           1.0 / 3.0 * C) +
      a_SM * B * 2 * alpha_ *
          square(baryon_effective_mass / electron_fraction) * A /
          rest_mass_density;
  const DataType rho_times_a_electron_second_log_density_deriv =
      a_electron_log_density_deriv * a_electron_log_density_deriv / a_electron +
      2.0 * a_electron / 9.0 * B * C * (1 - B * C);
  // One T already factored out
  const DataType degenerate =
      degenerate_temperature_dependent +
      temperature * (rho_times_a_SM_second_log_density_deriv +
                     rho_times_a_electron_second_log_density_deriv);

  return radiation + temperature * ideal * degenerate / (ideal + degenerate);
  ;
}

template <class ColdEos>
template <class DataType>
DataType AnalyticalThermal<ColdEos>::composition_dependent_internal_energy(
    const DataType& rest_mass_density,
    const DataType& electron_fraction) const {
  // Box 1
  const DataType equilibrium_electron_fraction =
      beta_equalibrium_proton_fraction(rest_mass_density);
  DataType result = 3.0 * K_ *
                    (cbrt(pow(electron_fraction, 4) * rest_mass_density) -
                     cbrt(pow(electron_fraction, 4) * rest_mass_density));
  result +=
      symmetry_energy_at_zero_temp(rest_mass_density) * 4 *
      (electron_fraction * (electron_fraction - 1.0) -
       equilibrium_electron_fraction * (equilibrium_electron_fraction - 1.0));
  return result;
}
template <class ColdEos>
template <class DataType>
DataType AnalyticalThermal<ColdEos>::composition_dependent_pressure(
    const DataType& rest_mass_density,
    const DataType& electron_fraction) const {
  const DataType equilibrium_electron_fraction =
      beta_equalibrium_proton_fraction(rest_mass_density);
  return K_ * rest_mass_density *
             (cbrt(pow(electron_fraction, 4) * rest_mass_density) -
              cbrt(pow(equilibrium_electron_fraction, 4) * rest_mass_density)) +
         symmetry_pressure_at_zero_temp(rest_mass_density) * 4 *
             (electron_fraction * (electron_fraction - 1.0) -
              equilibrium_electron_fraction *
                  (equilibrium_electron_fraction - 1.0));
}

template <typename ColdEos>
template <typename DataType>
DataType AnalyticalThermal<ColdEos>::dirac_effective_mass(
    const DataType& rest_mass_density) const {
  return invsqrt(1 + pow((rest_mass_density / n0_), -alpha_ * 2));
}

template <typename ColdEos>
template <typename DataType>
DataType AnalyticalThermal<ColdEos>::symmetry_energy_at_zero_temp(
    const DataType& rest_mass_density) const {
  // Erratum, Eq. (17)
  double common_factor = 0.6 * (pow(.25, 1.0 / 3.0) - 1.0);
  DataType kinetic_symmetry_component =
      common_factor * baryonic_fermi_internal_energy(rest_mass_density);
  double kinetic_symmetry_component_at_saturation =
      common_factor * baryonic_fermi_internal_energy(saturation_density_);
  // Eq. (14)
  return eta_ * kinetic_symmetry_component +
         (S0_ - eta_ * kinetic_symmetry_component_at_saturation) *
             pow(rest_mass_density / saturation_density_, gamma_);
}
template <class ColdEos>
template <class DataType>
DataType AnalyticalThermal<ColdEos>::beta_equalibrium_proton_fraction(
    const DataType& rest_mass_density) const {
  // Eq. (20) recast with z = (1 - 2 Y_{p,beta}) can be solved
  // via Cardano's formula.  I don't know if this is fast
  // but it's certainly simpler than unpacking and rootfinding
  DataType A = 64 / (3 * square(M_PI) * rest_mass_density) *
               cube(symmetry_energy_at_zero_temp(rest_mass_density) /
                    hbar_over_baryon_mass_to_four_thirds_);
  // z^3 + pz + q = 0; p = 1/(2A), q = -1/(2A)
  DataType p = 1 / (2 * A);
  // q = -p
  DataType sqrtDelta = p * sqrt(0.25 + p / 27.0);
  return cbrt(p * 0.5 + sqrtDelta) + cbrt(p * 0.5 - sqrtDelta);
}
template <class ColdEos>
template <typename DataType, typename MassType>
DataType AnalyticalThermal<ColdEos>::a_degeneracy(
    const DataType& rest_mass_density, const DataType& particle_fraction,
    const MassType& mass) const {
  const DataType kinetic_common_factor =
      pow(3 * square(M_PI) * particle_fraction * rest_mass_density, 2.0 / 3.0) *
      square(hbar_over_baryon_mass_to_four_thirds_);
  return square(M_PI) / 2 * sqrt(kinetic_common_factor + square(mass)) /
         kinetic_common_factor;
}
template <class ColdEos>
template <class DataType, class MassType>
DataType AnalyticalThermal<ColdEos>::a_degeneracy_log_density_derivative(
    const DataType& rest_mass_density, const DataType& particle_fraction,
    const MassType& mass) const {
  DataType log_mass_derivative_log_density =
      -alpha_ * (1 - square(mass / particle_fraction));
  DataType kinetic_common_factor =
      pow(3.0 * square(M_PI) * particle_fraction * rest_mass_density,
          2.0 / 3.0) *
      square(hbar_over_baryon_mass_to_four_thirds_ / mass);

  return 2.0 * a_degeneracy(rest_mass_density, particle_fraction, mass) / 3 *
         (1.0 -
          .5 * (1.0 / (1.0 + kinetic_common_factor) *
                (kinetic_common_factor + log_mass_derivative_log_density)));
}
template <class ColdEos>
template <class DataType>
DataType AnalyticalThermal<ColdEos>::symmetry_pressure_at_zero_temp(
    const DataType& rest_mass_density) const {
  double common_factor = 0.6 * (pow(.25, 1.0 / 3.0) - 1.0);
  DataType kinetic_symmetry_component =
      common_factor * baryonic_fermi_internal_energy(rest_mass_density);
  double kinetic_symmetry_component_at_saturation =
      common_factor * baryonic_fermi_internal_energy(saturation_density_);
  return rest_mass_density *
         (2.0 * eta_ / 3.0 * kinetic_symmetry_component +
          (S0_ - eta_ * kinetic_symmetry_component_at_saturation) *
              pow(rest_mass_density / saturation_density_, gamma_));
}
template <class ColdEos>
template <class DataType>
DataType
AnalyticalThermal<ColdEos>::symmetry_pressure_density_derivative_at_zero_temp(
    const DataType& rest_mass_density) const {
  double common_factor = 0.6 * (pow(.25, 1.0 / 3.0) - 1.0);
  DataType kinetic_symmetry_component =
      common_factor * baryonic_fermi_internal_energy(rest_mass_density);
  double kinetic_symmetry_component_at_saturation =
      common_factor * baryonic_fermi_internal_energy(saturation_density_);
  return rest_mass_density *
         (2.0 * eta_ * kinetic_symmetry_component +
          (S0_ - eta_ * kinetic_symmetry_component_at_saturation) *
              (gamma_ * (1 + gamma_)) *
              pow(rest_mass_density / saturation_density_, gamma_));
}

template <typename ColdEquationOfState>
template <class DataType>
Scalar<DataType> AnalyticalThermal<ColdEquationOfState>::
    pressure_from_density_and_temperature_impl(
        const Scalar<DataType>& rest_mass_density,
        const Scalar<DataType>& temperature,
        const Scalar<DataType>& electron_fraction) const {
  return Scalar<DataType>{
      get(cold_eos_.pressure_from_density(rest_mass_density)) +
      composition_dependent_pressure(get(rest_mass_density),
                                     get(electron_fraction)) +
      thermal_pressure(get(rest_mass_density), get(temperature),
                       get(electron_fraction))};
}

template <typename ColdEquationOfState>
template <class DataType>
Scalar<DataType>
AnalyticalThermal<ColdEquationOfState>::pressure_from_density_and_energy_impl(
    const Scalar<DataType>& rest_mass_density,
    const Scalar<DataType>& specific_internal_energy,
    const Scalar<DataType>& electron_fraction) const {
  // Have to rootfind for temperature
  const Scalar<DataType> temperature = temperature_from_density_and_energy(
      rest_mass_density, specific_internal_energy, electron_fraction);
  return pressure_from_density_and_temperature(rest_mass_density, temperature,
                                               electron_fraction);
}

template <typename ColdEquationOfState>
template <class DataType>
Scalar<DataType> AnalyticalThermal<ColdEquationOfState>::
    specific_internal_energy_from_density_and_temperature_impl(
        const Scalar<DataType>& rest_mass_density,
        const Scalar<DataType>& temperature,
        const Scalar<DataType>& electron_fraction) const {
  return Scalar<DataType>{
      get(cold_eos_.specific_internal_energy_from_density(rest_mass_density)) +
      composition_dependent_internal_energy(get(rest_mass_density),
                                            get(electron_fraction)) +
      thermal_internal_energy(get(rest_mass_density), get(temperature),
                              get(electron_fraction))};
}

template <typename ColdEquationOfState>
template <class DataType>
Scalar<DataType> AnalyticalThermal<ColdEquationOfState>::
    temperature_from_density_and_energy_impl(
        const Scalar<DataType>& rest_mass_density,
        const Scalar<DataType>& specific_internal_energy,
        const Scalar<DataType>& electron_fraction) const {
  const auto cold_energy =
      get(cold_eos_.specific_internal_energy_from_density(rest_mass_density));
  const auto composition_energy = composition_dependent_internal_energy(
      get(rest_mass_density), get(electron_fraction));
  const DataType thermal_energy_needed =
      get(specific_internal_energy) - cold_energy - composition_energy;
  // TODO: Do the right thing here pointwise
  const double fs = 11.0 / 4.0;
  const DataType radiation_prefactor =
      4 * stefan_boltzmann_sigma_ * fs / get(rest_mass_density);
  double ideal_prefactor = 1.5;
  const DataType degenerate_prefactor =
      (a_degeneracy(get(rest_mass_density),
                    make_with_value<DataType>(get(rest_mass_density), 0.5),
                    dirac_effective_mass(get(rest_mass_density))) +
       a_degeneracy(get(rest_mass_density), get(electron_fraction),
                    electron_mass_over_baryon_mass_));
  // Coefficients on rootfindng of a polynomial
  // E_thermal(T) = A_5 T^5 + A_4 T^4  + A_2 T^2 + A_1 T + A_0
  const DataType A0 = -thermal_energy_needed * ideal_prefactor;
  const DataType A1 = -thermal_energy_needed * degenerate_prefactor;
  const DataType A2 = ideal_prefactor * degenerate_prefactor;
  DataType A3;
  const DataType A4 = radiation_prefactor * ideal_prefactor;
  const DataType A5 = radiation_prefactor * degenerate_prefactor;
  // Most of the time the temperature will be small ( T < 1)
  double loose_upper_bound = LIKELY(max(thermal_energy_needed)) < 1.0e2
                                 ? 200.0
                                 : std::numeric_limits<double>::max();

  if constexpr (std::is_same_v<DataType, DataVector>) {
    A3 = make_with_value<DataVector>(get(rest_mass_density), 0.0);
  } else {
    A3 = 0.0;
  }
  return Scalar<DataType>{root_polynomial_positive(
      std::vector<DataType>{{A0, A1, A2, A3, A4, A5}}, 0.0, loose_upper_bound)};
}

template <typename ColdEquationOfState>
template <class DataType>
Scalar<DataType> AnalyticalThermal<ColdEquationOfState>::
    sound_speed_squared_from_density_and_temperature_impl(
        const Scalar<DataType>& rest_mass_density,
        const Scalar<DataType>& temperature,
        const Scalar<DataType>& electron_fraction) const {
  const DataType pressure_derivative_cold =
      get(cold_eos_.chi_from_density(rest_mass_density));
  const DataType equilibrium_electron_fraction =
      beta_equalibrium_proton_fraction(get(rest_mass_density));
  // Derived by hand, most likely to be wrong
  const DataType pressure_derivative_composition =
      (K_ * 4.0 / 3.0) *
          (cbrt(pow(get(electron_fraction), 4) * get(rest_mass_density)) -
           cbrt(pow(equilibrium_electron_fraction, 4) *
                get(rest_mass_density))) +
      symmetry_pressure_density_derivative_at_zero_temp(
          get(rest_mass_density)) *
          4.0 *
          (get(electron_fraction) * (get(electron_fraction) - 1.0) -
           equilibrium_electron_fraction *
               (equilibrium_electron_fraction - 1.0));

  const DataType thermal_pressure_derivative =
      thermal_pressure_density_derivative(
          get(rest_mass_density), get(temperature), get(electron_fraction));
  return Scalar<DataType>{
      1 /
      (1.0 +
       get(pressure_from_density_and_temperature(rest_mass_density, temperature,
                                                 electron_fraction)) /
           get(rest_mass_density) +
       get(specific_internal_energy_from_density_and_temperature(
           rest_mass_density, temperature, electron_fraction))) *
      (pressure_derivative_cold + pressure_derivative_composition +
       thermal_pressure_derivative)};
}
}  // namespace EquationsOfState

template class EquationsOfState::AnalyticalThermal<
    EquationsOfState::PolytropicFluid<true>>;
template class EquationsOfState::AnalyticalThermal<
    EquationsOfState::Enthalpy<EquationsOfState::Spectral>>;
