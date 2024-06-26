// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <boost/preprocessor/arithmetic/dec.hpp>
#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/control/expr_iif.hpp>
#include <boost/preprocessor/list/adt.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/tuple/to_list.hpp>
#include <limits>
#include <pup.h>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"  // IWYU pragma: keep
#include "PointwiseFunctions/Hydro/Units.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
/// \endcond

// IWYU pragma: no_forward_declare Tensor

namespace EquationsOfState {
/*!
 * \ingroup EquationsOfStateGroup
 * \brief Equation of state for a polytropic fluid
 *
 * A Cold Eos with an additional prescription to compute the electron fraction
 * as a function of density.  By doing this we avoid needing to specialize each
 * 1-D EoS to have reasonable electron fraction when constructing initla data.
 */
template <typename ColdEos, typename ElectronFractionComputer>
class Equilibrated1D : public EquationOfState<ColdEos::is_relativistic, 1> {
 public:
  static constexpr size_t thermodynamic_dim = 1;
  static constexpr bool is_relativistic = ColdEos::is_relativistic;

  struct UnderlyingEos {
    using type = ColdEos;
    static constexpr Options::String help = {
        "The 1-D EoS to augment with an electron fraction"};
  };
  struct ElectronFraction3DEos {
    using type = ElectronFractionComputer;
    static constexpr Options::String help = {
        "The 3-D Eos to take the electron fraction from in equilibrium"};
  };

  static constexpr Options::String help = {
      "A 1-d equation of state representing a 3D EoS in beta equilibrium at "
      "zero temperature"};

  using options = tmpl::list<UnderlyingEos, ElectronFraction3DEos>;

  Equilibrated1D() = default;
  Equilibrated1D(const Equilibrated1D&) = default;
  Equilibrated1D& operator=(const Equilibrated1D&) = default;
  Equilibrated1D(Equilibrated1D&&) = default;
  Equilibrated1D& operator=(Equilibrated1D&&) = default;
  ~Equilibrated1D() override = default;

  Equilibrated1D(const ColdEos& underlying_eos,
                 const ElectronFractionComputer& electron_fraction_computer);

  std::unique_ptr<EquationOfState<ColdEos::is_relativistic, 1>> get_clone()
      const override;

  std::unique_ptr<EquationOfState<ColdEos::is_relativistic, 3>>
  promote_to_3d_eos() const override;

  std::unique_ptr<EquationOfState<ColdEos::is_relativistic, 2>>
  promote_to_2d_eos() const override;

  bool is_equal(
      const EquationOfState<ColdEos::is_relativistic, 1>& rhs) const override;

  bool operator==(
      const Equilibrated1D<ColdEos, ElectronFractionComputer>& rhs) const;

  bool operator!=(
      const Equilibrated1D<ColdEos, ElectronFractionComputer>& rhs) const;

  Scalar<DataVector> equilibrium_electron_fraction_from_density_temperature(
      const Scalar<DataVector>& rest_mass_density,
      const Scalar<DataVector>& temperature) const override;

  Scalar<double> equilibrium_electron_fraction_from_density_temperature(
      const Scalar<double>& rest_mass_density,
      const Scalar<double>& temperature) const override;

  EQUATION_OF_STATE_FORWARD_DECLARE_MEMBERS(Equilibrated1D, 1)

  WRAPPED_PUPable_decl_base_template(  // NOLINT
      SINGLE_ARG(EquationOfState<is_relativistic, 1>), Equilibrated1D);

  /// The lower bound of the rest mass density that is valid for this EOS
  double rest_mass_density_lower_bound() const override;

  /// The upper bound of the rest mass density that is valid for this EOS
  double rest_mass_density_upper_bound() const override;

  /// The lower bound of the specific internal energy that is valid for this EOS
  /// at the given rest mass density \f$\rho\f$
  double specific_internal_energy_lower_bound(
      const double rest_mass_density) const override;

  /// The upper bound of the specific internal energy that is valid for this EOS
  /// at the given rest mass density \f$\rho\f$
  double specific_internal_energy_upper_bound(
      const double rest_mass_density) const override;

  /// The lower bound of the specific enthalpy that is valid for this EOS
  double specific_enthalpy_lower_bound() const override;

  /// The vacuum baryon mass for this EoS
  double baryon_mass() const override;

 private:
  EQUATION_OF_STATE_FORWARD_DECLARE_MEMBER_IMPLS(1)
  ColdEos underlying_eos_{};
  ElectronFractionComputer electron_fraction_computer_{};
};

/// \cond
template <typename ColdEos, typename ElectronFractionComputer>
PUP::able::PUP_ID EquationsOfState::Equilibrated1D<
    ColdEos, ElectronFractionComputer>::my_PUP_ID = 0;
/// \endcond
}  // namespace EquationsOfState
