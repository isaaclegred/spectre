# Distributed under the MIT License.
# See LICENSE.txt for details.

from PolytropicFluid import polytropic_pressure_from_density
from PolytropicFluid import polytropic_specific_internal_energy_from_density


class PolytropicFluid:
    def __init__(self, polytropic_constant, polytropic_exponent):
        self.polytropic_constant_ = polytropic_constant
        self.polytropic_exponent_ = polytropic_exponent

    def pressure_from_density(self, rho):
        return self.polytropic_constant_ * rho**(polytropic_exponent_)

    def internal_energy_from_density(self, rho):
        return pressure_from_density(rho) / (rho *
                                             (self.polytropic_constant_ - 1.0))

    def chi_from_density(self, rho):
        return polytropic_constant_ * polytropic_exponent_ * rho**(
            polytropic_exponent_ - 1.0)


class AnalyticalThermal:
    def __init__(self, S0, L, gamma, n0, alpha, cold_eos):
        self.S0_ = S0
        self.L_ = L
        self.gamma_ = gamma
        self.n0_ = n0
        self.alpha_ = alpha
        self.cold_eos_ = cold_eos


def analytical_thermal_polytrope_pressure_from_density_and_energy(
    rest_mass_density, specific_internal_energy, polytropic_constant,
    polytropic_exponent, thermal_adiabatic_index):
    p_c = polytropic_pressure_from_density(rest_mass_density,
                                           polytropic_constant,
                                           polytropic_exponent)
    eps_c = polytropic_specific_internal_energy_from_density(
        rest_mass_density, polytropic_constant, polytropic_exponent)
    return p_c + rest_mass_density * (specific_internal_energy -
                                      eps_c) * (thermal_adiabatic_index - 1.0)


def hybrid_polytrope_rel_pressure_from_density_and_enthalpy(
    rest_mass_density, specific_enthalpy, polytropic_constant,
    polytropic_exponent, thermal_adiabatic_index):
    p_c = polytropic_pressure_from_density(rest_mass_density,
                                           polytropic_constant,
                                           polytropic_exponent)
    eps_c = polytropic_specific_internal_energy_from_density(
        rest_mass_density, polytropic_constant, polytropic_exponent)
    return p_c / thermal_adiabatic_index + (rest_mass_density *
                                            (specific_enthalpy - 1.0 - eps_c) *
                                            (thermal_adiabatic_index - 1.0) /
                                            thermal_adiabatic_index)


def hybrid_polytrope_temperature_from_density_and_energy(
    rest_mass_density, specific_internal_energy, polytropic_constant,
    polytropic_exponent, thermal_adiabatic_index):
    eps_c = polytropic_specific_internal_energy_from_density(
        rest_mass_density, polytropic_constant, polytropic_exponent)
    return (thermal_adiabatic_index - 1.0) * (specific_internal_energy - eps_c)


def hybrid_polytrope_newt_pressure_from_density_and_enthalpy(
    rest_mass_density, specific_enthalpy, polytropic_constant,
    polytropic_exponent, thermal_adiabatic_index):
    p_c = polytropic_pressure_from_density(rest_mass_density,
                                           polytropic_constant,
                                           polytropic_exponent)
    eps_c = polytropic_specific_internal_energy_from_density(
        rest_mass_density, polytropic_constant, polytropic_exponent)
    return p_c / thermal_adiabatic_index + (rest_mass_density *
                                            (specific_enthalpy - eps_c) *
                                            (thermal_adiabatic_index - 1.0) /
                                            thermal_adiabatic_index)


def hybrid_polytrope_rel_specific_enthalpy_from_density_and_energy(
    rest_mass_density, specific_internal_energy, polytropic_constant,
    polytropic_exponent, thermal_adiabatic_index):
    p_c = polytropic_pressure_from_density(rest_mass_density,
                                           polytropic_constant,
                                           polytropic_exponent)
    eps_c = polytropic_specific_internal_energy_from_density(
        rest_mass_density, polytropic_constant, polytropic_exponent)
    return 1.0 + eps_c + p_c / rest_mass_density + thermal_adiabatic_index * (
        specific_internal_energy - eps_c)


def hybrid_polytrope_newt_specific_enthalpy_from_density_and_energy(
    rest_mass_density, specific_internal_energy, polytropic_constant,
    polytropic_exponent, thermal_adiabatic_index):
    p_c = polytropic_pressure_from_density(rest_mass_density,
                                           polytropic_constant,
                                           polytropic_exponent)
    eps_c = polytropic_specific_internal_energy_from_density(
        rest_mass_density, polytropic_constant, polytropic_exponent)
    return eps_c + p_c / rest_mass_density + thermal_adiabatic_index * (
        specific_internal_energy - eps_c)


def hybrid_polytrope_specific_internal_energy_from_density_and_pressure(
    rest_mass_density, pressure, polytropic_constant, polytropic_exponent,
    thermal_adiabatic_index):
    p_c = polytropic_pressure_from_density(rest_mass_density,
                                           polytropic_constant,
                                           polytropic_exponent)
    eps_c = polytropic_specific_internal_energy_from_density(
        rest_mass_density, polytropic_constant, polytropic_exponent)
    return eps_c + (pressure - p_c) / (thermal_adiabatic_index -
                                       1.0) / rest_mass_density


def hybrid_polytrope_chi_from_density_and_energy(rest_mass_density,
                                                 specific_internal_energy,
                                                 polytropic_constant,
                                                 polytropic_exponent,
                                                 thermal_adiabatic_index):
    p_c = polytropic_pressure_from_density(rest_mass_density,
                                           polytropic_constant,
                                           polytropic_exponent)
    eps_c = polytropic_specific_internal_energy_from_density(
        rest_mass_density, polytropic_constant, polytropic_exponent)
    return polytropic_exponent * p_c / rest_mass_density + (
        specific_internal_energy - eps_c -
        p_c / rest_mass_density) * (thermal_adiabatic_index - 1.0)


def hybrid_polytrope_kappa_times_p_over_rho_squared_from_density_and_energy(
    rest_mass_density, specific_internal_energy, polytropic_constant,
    polytropic_exponent, thermal_adiabatic_index):
    p_c = polytropic_pressure_from_density(rest_mass_density,
                                           polytropic_constant,
                                           polytropic_exponent)
    eps_c = polytropic_specific_internal_energy_from_density(
        rest_mass_density, polytropic_constant, polytropic_exponent)
    return (thermal_adiabatic_index - 1.0) * p_c / rest_mass_density + (
        specific_internal_energy - eps_c) * (thermal_adiabatic_index - 1.0)**2
