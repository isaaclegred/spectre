// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/Hydro/EquationsOfState/Equilibrated1D.hpp"

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <cmath>
#include <limits>
#include <memory>

#include "DataStructures/DataVector.hpp"  // IWYU pragma: keep
#include "DataStructures/Tensor/Tensor.hpp"

#include "PointwiseFunctions/Hydro/EquationsOfState/AnalyticalThermal.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/Barotropic2D.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/Barotropic3D.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/PolytropicFluid.hpp"
#include "Utilities/MakeWithValue.hpp"

// IWYU pragma: no_forward_declare Tensor

namespace EquationsOfState {
template <typename ColdEos, typename ElectronFractionComputer>
Equilibrated1D<ColdEos, ElectronFractionComputer>::Equilibrated1D(
    const ColdEos& underlying_eos,
    const ElectronFractionComputer& electron_fraction_computer)
    : underlying_eos_(underlying_eos),
      electron_fraction_computer_(electron_fraction_computer) {}

EQUATION_OF_STATE_MEMBER_DEFINITIONS(
    SINGLE_ARG(template <typename ColdEos, typename ElectronFractionComputer>),
    SINGLE_ARG(Equilibrated1D<ColdEos, ElectronFractionComputer>), double, 1)
EQUATION_OF_STATE_MEMBER_DEFINITIONS(
    SINGLE_ARG(template <typename ColdEos, typename ElectronFractionComputer>),
    SINGLE_ARG(Equilibrated1D<ColdEos, ElectronFractionComputer>), DataVector,
    1)

template <typename ColdEos, typename ElectronFractionComputer>
bool Equilibrated1D<ColdEos, ElectronFractionComputer>::operator==(
    const Equilibrated1D<ColdEos, ElectronFractionComputer>& rhs) const {
  return underlying_eos_ == rhs.underlying_eos_ and
         electron_fraction_computer_ == rhs.electron_fraction_computer_;
}

template <typename ColdEos, typename ElectronFractionComputer>
bool Equilibrated1D<ColdEos, ElectronFractionComputer>::operator!=(
    const Equilibrated1D<ColdEos, ElectronFractionComputer>& rhs) const {
  return not(*this == rhs);
}

template <typename ColdEos, typename ElectronFractionComputer>
bool Equilibrated1D<ColdEos, ElectronFractionComputer>::is_equal(
    const EquationOfState<ColdEos::is_relativistic, 1>& rhs) const {
  const auto& derived_ptr = dynamic_cast<
      const Equilibrated1D<ColdEos, ElectronFractionComputer>* const>(&rhs);
  return derived_ptr != nullptr and *derived_ptr == *this;
}

template <typename ColdEos, typename ElectronFractionComputer>
std::unique_ptr<EquationOfState<ColdEos::is_relativistic, 1>>
Equilibrated1D<ColdEos, ElectronFractionComputer>::get_clone() const {
  auto clone =
      std::make_unique<Equilibrated1D<ColdEos, ElectronFractionComputer>>(
          *this);
  return std::unique_ptr<EquationOfState<is_relativistic, 1>>(std::move(clone));
}

template <typename ColdEos, typename ElectronFractionComputer>
std::unique_ptr<EquationOfState<ColdEos::is_relativistic, 3>>
Equilibrated1D<ColdEos, ElectronFractionComputer>::promote_to_3d_eos() const {
  return std::make_unique<
      Barotropic3D<Equilibrated1D<ColdEos, ElectronFractionComputer>>>(*this);
}

template <typename ColdEos, typename ElectronFractionComputer>
std::unique_ptr<EquationOfState<ColdEos::is_relativistic, 2>>
Equilibrated1D<ColdEos, ElectronFractionComputer>::promote_to_2d_eos() const {
  return std::make_unique<
      Barotropic2D<Equilibrated1D<ColdEos, ElectronFractionComputer>>>(*this);
}
template <typename ColdEos, typename ElectronFractionComputer>
Equilibrated1D<ColdEos, ElectronFractionComputer>::Equilibrated1D(
    CkMigrateMessage* msg)
    : EquationOfState<is_relativistic, 1>(msg) {}

template <typename ColdEos, typename ElectronFractionComputer>
void Equilibrated1D<ColdEos, ElectronFractionComputer>::pup(PUP::er& p) {
  EquationOfState<is_relativistic, 1>::pup(p);
  p | underlying_eos_;
  p | electron_fraction_computer_;
}

template <typename ColdEos, typename ElectronFractionComputer>
template <class DataType>
Scalar<DataType>
Equilibrated1D<ColdEos, ElectronFractionComputer>::pressure_from_density_impl(
    const Scalar<DataType>& rest_mass_density) const {
  return underlying_eos_.pressure_from_density(rest_mass_density);
}

template <typename ColdEos, typename ElectronFractionComputer>
template <class DataType>
Scalar<DataType> Equilibrated1D<ColdEos, ElectronFractionComputer>::
    rest_mass_density_from_enthalpy_impl(
        const Scalar<DataType>& specific_enthalpy) const {
  return underlying_eos_.rest_mass_density_from_enthalpy(specific_enthalpy);
}

template <typename ColdEos, typename ElectronFractionComputer>
template <class DataType>
Scalar<DataType> Equilibrated1D<ColdEos, ElectronFractionComputer>::
    specific_internal_energy_from_density_impl(
        const Scalar<DataType>& rest_mass_density) const {
  return underlying_eos_.specific_internal_energy_from_density(
      rest_mass_density);
}

template <typename ColdEos, typename ElectronFractionComputer>
template <class DataType>
Scalar<DataType>
Equilibrated1D<ColdEos, ElectronFractionComputer>::chi_from_density_impl(
    const Scalar<DataType>& rest_mass_density) const {
  return underlying_eos_.chi_from_density(rest_mass_density);
}

template <typename ColdEos, typename ElectronFractionComputer>
template <class DataType>
Scalar<DataType> Equilibrated1D<ColdEos, ElectronFractionComputer>::
    kappa_times_p_over_rho_squared_from_density_impl(
        const Scalar<DataType>& rest_mass_density) const {
  return make_with_value<Scalar<DataType>>(get(rest_mass_density), 0.0);
}

template <typename ColdEos, typename ElectronFractionComputer>
double Equilibrated1D<
    ColdEos, ElectronFractionComputer>::rest_mass_density_upper_bound() const {
  return underlying_eos_.rest_mass_density_upper_bound();
}

template <typename ColdEos, typename ElectronFractionComputer>
double Equilibrated1D<
    ColdEos, ElectronFractionComputer>::rest_mass_density_lower_bound() const {
  return underlying_eos_.rest_mass_density_lower_bound();
}
template <typename ColdEos, typename ElectronFractionComputer>
double Equilibrated1D<ColdEos, ElectronFractionComputer>::
    specific_internal_energy_lower_bound(const double rest_mass_density) const {
  return underlying_eos_.specific_internal_energy_lower_bound(
      rest_mass_density);
}

template <typename ColdEos, typename ElectronFractionComputer>
double Equilibrated1D<ColdEos, ElectronFractionComputer>::
    specific_internal_energy_upper_bound(const double rest_mass_density) const {
  return underlying_eos_.specific_internal_energy_upper_bound(
      rest_mass_density);
}

template <typename ColdEos, typename ElectronFractionComputer>
double Equilibrated1D<
    ColdEos, ElectronFractionComputer>::specific_enthalpy_lower_bound() const {
  return underlying_eos_.specific_enthalpy_lower_bound();
}

template <typename ColdEos, typename ElectronFractionComputer>
double Equilibrated1D<ColdEos, ElectronFractionComputer>::baryon_mass() const {
  return underlying_eos_.baryon_mass();
}

template <typename ColdEos, typename ElectronFractionComputer>
Scalar<DataVector> Equilibrated1D<ColdEos, ElectronFractionComputer>::
    equilibrium_electron_fraction_from_density_temperature(
        const Scalar<DataVector>& rest_mass_density,
        const Scalar<DataVector>& temperature) const {
  return electron_fraction_computer_
      .equilibrium_electron_fraction_from_density_temperature(rest_mass_density,
                                                              temperature);
}

template <typename ColdEos, typename ElectronFractionComputer>

Scalar<double> Equilibrated1D<ColdEos, ElectronFractionComputer>::
    equilibrium_electron_fraction_from_density_temperature(
        const Scalar<double>& rest_mass_density,
        const Scalar<double>& temperature) const {
  return electron_fraction_computer_
      .equilibrium_electron_fraction_from_density_temperature(rest_mass_density,
                                                              temperature);
}
}  // namespace EquationsOfState

template class EquationsOfState::Equilibrated1D<
    EquationsOfState::PolytropicFluid<true>,
    EquationsOfState::AnalyticalThermal<
        EquationsOfState::PolytropicFluid<true>>>;
// template class EquationsOfState::Equilibrated1D<false>;
