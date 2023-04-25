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
#include <cstddef>
#include <limits>
#include <pup.h>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/Options.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
/// \endcond

namespace EquationsOfState {
/*!
 * \ingroup EquationsOfStateGroup
 *
 * \brief Analytical extension of a cold EoS to arbitrary temperature and
 * proton fraction using the expressions in Raithel, Özel, and and Psaltis
 * DOI : 10.3847/1538-4357/ab08ea
 *
 * The analytical extension of the EoS includes both thermal and
 * out-of-beta-equalibrium corrections to the EoS.  The extension proceeds
 * by expressing the energy per baryon \f$E/A = m_b + \epsilon m_b\f$ as a sum
 of a
 * cold part (at temperature \f$T=0\f$ and at proton fraction in beta
 * equalibrium \f$Y_p = Y_{p,\beta})\f$, a component
 * dependent on electron (or proton) fraction and
 * baryon density, and finally a component dependent on temperature,
 * density, and composition.
 *
 * Box I gives the total energy per particle
 * \f[
 * E/A(n, Y_p, T) = E/A(n, Y_{p, \beta}, T=0) +
 * 3K(Y_p^{4/3} - Y_{p, \beta}^4/3)n^{1/3} +
 * E_{\rm sym}(n, T=0)\left[(1-2Y_p)^2 - (1-2Y_{p,\beta})^2 +
 * E_{\rm thm}(n, Y_p, T)\right]
 * \f]
 *
 * where  \f$\n\f$ is the number density of baryons, and
 * \f[
 * E_{\rm thm}(n, Y_p, T) =
 *  4 \sigma f_s T^4/(c n^4)  +
 *  \left{[(3/2) k_B T]^{-1}
 *  [a(n, M^*_{\rm SM})T^2 + a(n, Y_p, m_e)Y_pT^2 ]^-1\right}^{-1}
 * \f]
 *  is the thermal part of the energy.
 *
 * The various components of the energy are defined as
 * \f$K \equiv (3\pi^2)^{1/3}(\hbar c / 4)\f$
 * \f[
 * E_{\rm sym}(m, T=0) = \eta E^{\rm kin}_{\rm sym}(n) +
 * [S_0-\eta E^{\rm kin}_{\rm sym}(n_sat) ]\left(\frac{n}{n_{\rm
 sat}}\right)^{\gamma}
 * \f]
 * \f[
 * E^{\rm kin}_{\rm sym} = \frac{3}{5}(2^{-2/3}-1 )E_f(n)
 * \f]
 * (Note this equation is corrected in an erratum 10.3847/1538-4357/ac0630)
 * \f[
 * \eta = \frac{5/9}\left[\frac{L - 3 S_0 \gamma}{(2^{-2/3}-1)(2/3-\gamma)
 * E_f(n_sat)}\right]
 * \f]
 * And \f$E_f(n)\f$ is the Fermi energy
 * \f[
 * E_f(n) = \frac{\hbar^2}{2m }(3\pi^2 n)^{2/3}
 * \f]
 * \f$S_0, L\f$ are the symmetry energy and the slope of the symmetry energy
 * at saturation density, respectively,
 * \f$\gamma\f$ is a constant dependent on the EoS
 * \f$m_e\f$ is the mass of the electron
 * \f$ \sigma, f_s\f$ are the Stefan-Boltzmann constant and the
 * See the erratum for a more explicit expression of \f$a\f$
 * \f[
 * a(n,Y_p, m) = \frac{\pi^2k_B^2}/2 \frac{\sqrt{(3\pi^2Y_pn)^{2/3}
 * (\hbar c^2) + m^2 }}{(3\pi^2Y_p n)^{2/3}(\hbar c)^2}
 * \f]
 * The dirac effective mass is
 * \f[
 * M^{*}_{\rm SM}\left{(mc^2)^{-1/2} +
 \left[mc^2\left(\frac{n}{n_0}\right)^{-\alpha}\right]^{-2}\right}^{-1/2}
 * \f]
 * and \f$\alpha\f$ is a parameter which controls the rate of fall off of the
 effective
 * mass.
 * The scheme additionally includes a strategy for computing $Y_{p, \beta}$.
 * We omit the expression here because we do not calculate using it,
 *  but we note it depends on
 * \f$L\f$, \f$S_0\f$ and $\gamma$, and satisfies
 * \f[
 * Y_{p, \beta} = (1-2Y_{p, \beta}) \frac{64}{3 pi^2 n}\left[\frac{E_{\rm sym
 (n, T=0)}}{\hbar x}\right]^3
 * \f]
 *
 * Equations are framed in the paper in arbitrary units though computations are
 given
 * in nuclear units \f$[S_0] = [L] = \mathrm{MeV}\f$, \f$[n] = \mathrm{fm}^{-3}
 \f$
 * we work in geometric units so
 * that \f$G=c=M_{\odot}=1\f$.  Therefore \f$ \hbar \sim 1 \times 10^{-76}
 M_{\odot}^2,\
 * $n_{\rm sat} \sim 8 \time 10^{56} M_{\odot}^{-3} \f$.
 * This unit system is not well adapted, so instead we use
 * \f$ \rho_{\rm sat} = m_{N} * n_[\rm sat}] \sim 4 \times 10^{-4}}  \f$,
 * Likewise, we don't want to consider the energy per particle, we want the
 * specific energy, $1 + \epsilon$, so we divide all energies by the baryon
 mass.  I.e.
 * every energy should be given in dimensionless units representing fractions of
 a baryon mass.
 * For example, \f$ k_B = 1\f$ so all temperatures are measured in baryon
 masses.

 */
template <typename ColdEquationOfState>
class AnalyticalThermal
    : public EquationOfState<ColdEquationOfState::is_relativistic, 3> {
 public:
  static constexpr size_t thermodynamic_dim = 3;
  static constexpr bool is_relativistic = ColdEquationOfState::is_relativistic;

  struct ColdEos {
    using type = ColdEquationOfState;
    static constexpr Options::String help = {"Cold equation of state"};
  };

  struct SymmetryEnergy {
    using type = double;
    static constexpr Options::String help = {
        "The symmetry energy of EoS evaluated at saturation density, S_0"};
  };
  struct SymmetryEnergySlope {
    using type = double;
    static constexpr Options::String help = {
        "The slope of the symmetry"
        "energy at saturation density, L. In units of baryon mass"};
  };
  struct SymmetryEnergyGamma {
    using type = double;
    static constexpr Options::String help = {
        "Phenomenological parameter"
        "gamma controlling the dependence of the potential symmetry energy"
        "away from saturation density."};
  };
  struct EffectiveMassn0 {
    using type = double;
    static constexpr Options::String help = {
        "Phenomenological"
        "parameter, n_0, controlling the onset of the effective mass drop off "
        "."};
  };
  struct EffectiveMassAlpha {
    using type = double;
    static constexpr Options::String help = {
        "Phenomenological"
        "parameter, alpha,  controlling the rate of the effective mass drop "
        "off ."};
  };

  static constexpr Options::String help = {
      "A hybrid equation of state combining a cold EOS with a simple thermal "
      "part.  The pressure is related to the rest mass density by "
      " p = p_cold(rho) + rho * (Gamma_th - 1) * (epsilon - "
      "epsilon_cold(rho)), where p is the pressure, rho is the rest mass "
      "density, epsilon is the specific internal energy, p_cold and "
      "epsilon_cold are the pressure and specific internal energy evaluated "
      "using the cold EOS and Gamma_th is the adiabatic index for the thermal "
      "part."};

  using options =
      tmpl::list<ColdEos, SymmetryEnergy, SymmetryEnergySlope,
                 SymmetryEnergyGamma, EffectiveMassn0, EffectiveMassAlpha>;

  AnalyticalThermal() = default;
  AnalyticalThermal(const AnalyticalThermal&) = default;
  AnalyticalThermal& operator=(const AnalyticalThermal&) = default;
  AnalyticalThermal(AnalyticalThermal&&) = default;
  AnalyticalThermal& operator=(AnalyticalThermal&&) = default;
  ~AnalyticalThermal() override = default;

  AnalyticalThermal(ColdEquationOfState cold_eos, double S0, double L,
                    double gamma, double n0, double alpha);

  EQUATION_OF_STATE_FORWARD_DECLARE_MEMBERS(AnalyticalThermal, 3)

  WRAPPED_PUPable_decl_base_template(  // NOLINT
      SINGLE_ARG(EquationOfState<is_relativistic, 3>), AnalyticalThermal);

  std::unique_ptr<EquationOfState<is_relativistic, 3>> get_clone()
      const override;

  bool operator==(const AnalyticalThermal<ColdEquationOfState>& rhs) const;

  bool operator!=(const AnalyticalThermal<ColdEquationOfState>& rhs) const;

  bool is_equal(const EquationOfState<is_relativistic, 3>& rhs) const override;

  /// The lower bound of the electron fraction that is valid for this EOS
  double electron_fraction_lower_bound() const override { return 0.0; }

  /// The upper bound of the electron fraction that is valid for this EOS
  double electron_fraction_upper_bound() const override { return 0.5; }

  /// The lower bound of the temperature that is valid for this EOS
  double temperature_lower_bound() const override { return 0.0; };

  /// The upper bound of the temperature that is valid for this EOS
  double temperature_upper_bound() const override {
    return std::numeric_limits<double>::max();
  };

  /// The lower bound of the rest mass density that is valid for this EOS
  double rest_mass_density_lower_bound() const override {
    return cold_eos_.rest_mass_density_lower_bound();
  }

  /// The upper bound of the rest mass density that is valid for this EOS
  double rest_mass_density_upper_bound() const override {
    return cold_eos_.rest_mass_density_upper_bound();
  }

  /// The lower bound of the specific internal energy that is valid for this EOS
  /// at the given rest mass density \f$\rho\f$
  double specific_internal_energy_lower_bound(
      const double rest_mass_density,
      const double /*electron_fraction*/) const override {
    return cold_eos_.specific_internal_energy_lower_bound(rest_mass_density);
  }

  /// The upper bound of the specific internal energy that is valid for this EOS
  /// at the given rest mass density \f$\rho\f$
  double specific_internal_energy_upper_bound(
      const double /*rest_mass_density*/,
      const double /*electron_fraction*/) const override {
    return std::numeric_limits<double>::max();
  }

  /// The lower bound of the specific enthalpy that is valid for this EOS
  double specific_enthalpy_lower_bound() const override {
    return cold_eos_.specific_enthalpy_lower_bound();
  }

 private:
  template <class DataType>
  DataType baryonic_fermi_internal_energy(
      const DataType& rest_mass_density) const;
  template <class DataType>
  DataType thermal_internal_energy(const DataType& rest_mass_density,
                                   const DataType& temperature,
                                   const DataType& electron_fraction) const;
  template <class DataType>
  DataType thermal_pressure(const DataType& rest_mass_density,
                            const DataType& temperature,
                            const DataType& electron_fraction) const;
  template <class DataType>
  DataType thermal_pressure_density_derivative(
      const DataType& rest_mass_density, const DataType& temperature,
      const DataType& electron_fraction) const;
  template <class DataType>
  DataType composition_dependent_internal_energy(
      const DataType& rest_mass_density,
      const DataType& electron_fraction) const;
  template <class DataType>
  DataType composition_dependent_pressure(
      const DataType& rest_mass_density,
      const DataType& electron_fraction) const;
  template <class DataType>
  DataType dirac_effective_mass(const DataType& rest_mass_density) const;
  template <class DataType>
  DataType symmetry_energy_at_zero_temp(
      const DataType& rest_mass_density) const;
  template <class DataType>
  DataType beta_equalibrium_proton_fraction(
      const DataType& rest_mass_density) const;
  template <class DataType, typename MassType>
  DataType a_degeneracy(const DataType& rest_mass_density,
                        const DataType& electron_fraction,
                        const MassType& mass) const;
  template <class DataType, typename MassType>
  DataType a_degeneracy_log_density_derivative(
      const DataType& rest_mass_density, const DataType& electron_fraction,
      const MassType& mass) const;
  double get_eta() const;

  template <class DataType>
  DataType symmetry_pressure_at_zero_temp(
      const DataType& rest_mass_density) const;

  template <class DataType>
  DataType symmetry_pressure_density_derivative_at_zero_temp(
      const DataType& rest_mass_density) const;

  EQUATION_OF_STATE_FORWARD_DECLARE_MEMBER_IMPLS(3)

  ColdEquationOfState cold_eos_;
  double S0_;
  double L_;
  double gamma_;
  double n0_;
  double alpha_;
  double eta_;
  // This variable just happens to be the way that
  // these microscopic quantities enter all of the
  // EoS expressions. It is not dimensionless in
  // geometrized units
  double hbar_over_baryon_mass_to_four_thirds_ = 1.5060;
  double saturation_density_ = 4.53e-4;
  double stefan_boltzmann_sigma_ =
      2.0 * pow(M_PI, 5) / 15 *
      pow(1 / hbar_over_baryon_mass_to_four_thirds_, 3);
  double electron_mass_over_baryon_mass_ = 5.44e-4;
  double K_ =
      cbrt(3 * square(M_PI)) * hbar_over_baryon_mass_to_four_thirds_ * 0.25;
};

/// \cond
template <typename ColdEquationOfState>
PUP::able::PUP_ID
    EquationsOfState::AnalyticalThermal<ColdEquationOfState>::my_PUP_ID = 0;
/// \endcond
}  // namespace EquationsOfState
