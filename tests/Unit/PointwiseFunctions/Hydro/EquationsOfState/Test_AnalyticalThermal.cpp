// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Framework/TestingFramework.hpp"

#include <limits>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestCreation.hpp"
#include "Helpers/PointwiseFunctions/Hydro/EquationsOfState/TestHelpers.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/Factory.hpp"
#include "PointwiseFunctions/Hydro/SpecificEnthalpy.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
namespace {
// void check_random_polytrope() {
//   register_derived_classes_with_charm<
//       EquationsOfState::EquationOfState<true, 3>>();
//   const double d_for_size = std::numeric_limits<double>::signaling_NaN();
//   const DataVector dv_for_size(5);
//   TestHelpers::EquationsOfState::check(
//       EquationsOfState::HybridEos<
//           EquationsOfState::PolytropicFluid<IsRelativistic>>{
//           EquationsOfState::PolytropicFluid<IsRelativistic>{100.0, 4.0
//           / 3.0}, 5.0 / 3.0},
//       "HybridEos", "hybrid_polytrope", d_for_size, 100.0, 4.0 / 3.0, 5.0
//       / 3.0);
//   TestHelpers::EquationsOfState::check(
//       EquationsOfState::HybridEos<
//           EquationsOfState::PolytropicFluid<IsRelativistic>>{
//           EquationsOfState::PolytropicFluid<IsRelativistic>{100.0, 4.0
//           / 3.0}, 5.0 / 3.0},
//       "HybridEos", "hybrid_polytrope", dv_for_size, 100.0, 4.0 / 3.0,
//       5.0 / 3.0);
// }

void check_exact_polytrope() {
  EquationsOfState::PolytropicFluid<true> cold_eos{100.0, 2.0};
  const Scalar<double> rho_c{1.0e-3};
  const auto p_c = cold_eos.pressure_from_density(rho_c);
  CHECK(get(p_c) == approx(1.0e-4));
  const auto eps_c = cold_eos.specific_internal_energy_from_density(rho_c);
  CHECK(get(eps_c) == approx(1.0e-1));
  const auto h_c = hydro::relativistic_specific_enthalpy(rho_c, eps_c, p_c);
  CHECK(get(h_c) == approx(1.0 + 1.0e-1 + 1.0e-1));
  const auto chi_c = cold_eos.chi_from_density(rho_c);
  CHECK(get(chi_c) == approx(2.0e-1));
  const auto p_c_kappa_c_over_rho_sq =
      cold_eos.kappa_times_p_over_rho_squared_from_density(rho_c);
  CHECK(get(p_c_kappa_c_over_rho_sq) == 0.0);
  const auto c_s_sq = (get(chi_c) + get(p_c_kappa_c_over_rho_sq)) / get(h_c);
  CHECK(c_s_sq == approx(1.0 / 6.0));
  EquationsOfState::AnalyticalThermal<EquationsOfState::PolytropicFluid<true>>
      eos{{100.0, 2.0}, 1.5, .1, .1, 1.0, .85};
  TestHelpers::EquationsOfState::test_get_clone(eos);

  EquationsOfState::AnalyticalThermal<EquationsOfState::PolytropicFluid<true>>
      other_eos{{100.0, 3.0}, 1.5, .1, .1, 1.0, .89};
  const auto other_type_eos =
      EquationsOfState::PolytropicFluid<true>{100.0, 2.0};
  const Scalar<double> rho{.001};
  CHECK(eos == eos);
  CHECK(eos != other_eos);
  CHECK(eos != other_type_eos);
  const Scalar<double> eps{.2};
  const Scalar<double> temp{.01};
  const Scalar<double> temp_times_two{.02};
  const Scalar<double> elec_frac{.05};
  const Scalar<double> small_rho { 1.0e-6 }

  const auto p = eos.pressure_from_density_and_energy(rho, eps, elec_frac);
  CHECK(get(p) == approx(1e-4));
  CHECK(get(p) >= get(p_c));
  const auto h = hydro::relativistic_specific_enthalpy(rho, eps, p);
  CHECK(get(h) == approx(1.0));
  CHECK(get(h) >= get(h_c));
  CHECK(get(eos.pressure_from_density_and_temperature(rho, temp, elec_frac)) >=
        get(p_c));
  CHECK(get(eos.pressure_from_density_and_temperature(rho, temp_times_two,
                                                      elec_frac)) /
            square(get(temp_times_two)) ==
        approx(1.5));

  CHECK(get(eos.pressure_from_density_and_temperature(rho, temp_times_two,
                                                      elec_frac)) /
            square(get(temp_times_two)) ==
        approx(1.5));
  CHECK(get(eos.specific_internal_energy_from_density_and_temperature(
            rho, temp, elec_frac)) >= get(eps_c));
  const auto temp_recovered =
      eos.temperature_from_density_and_energy(rho, eps, elec_frac);
  CHECK(get(temp_recovered) == approx(.01));
  // Check consistency of eps(rho, T(rho, eps', Ye), Ye) = eps'
  CHECK(get(eos.specific_internal_energy_from_density_and_temperature(
            rho, temp_recovered, elec_frac)) == approx(get(eps)));
  const auto speed_of_sound =
      eos.sound_speed_squared_from_density_and_temperature(rho, eps, elec_frac);
  CHECK(get(speed_of_sound) == approx(.33));
}

void check_bounds() {
  const auto cold_eos = EquationsOfState::PolytropicFluid<true>{100.0, 1.5};
  const EquationsOfState::AnalyticalThermal<
      EquationsOfState::PolytropicFluid<true>>
      eos{cold_eos, 1.5, .1, .1, 1.0, .89};
  double electron_fraction = .05;
  double rest_mass_density = .005;
  CHECK(0.0 == eos.rest_mass_density_lower_bound());
  CHECK(0.0 == eos.temperature_lower_bound());

  CHECK(0.0 == eos.specific_internal_energy_lower_bound(rest_mass_density,
                                                        electron_fraction));
  CHECK(0.0 == eos.electron_fraction_lower_bound());
  CHECK(0.5 == eos.electron_fraction_upper_bound());

  const double max_double = std::numeric_limits<double>::max();
  CHECK(max_double == eos.rest_mass_density_upper_bound());
  CHECK(max_double == eos.specific_internal_energy_upper_bound(
                          rest_mass_density, electron_fraction));
  CHECK(max_double == eos.temperature_upper_bound());
}

}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.EquationsOfState.AnalyticalThermal",
                  "[Unit][EquationsOfState]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/Hydro/EquationsOfState/"};
  // check_random_polytrope<true>();
  check_exact_polytrope();
  check_bounds();
}
