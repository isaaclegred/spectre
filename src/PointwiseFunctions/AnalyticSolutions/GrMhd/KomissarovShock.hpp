// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <limits>
#include <string>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GrMhd/Solutions.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/IdealFluid.hpp"  // IWYU pragma: keep
#include "PointwiseFunctions/Hydro/TagsDeclarations.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// IWYU pragma: no_include <pup.h>

/// \cond
namespace PUP {
class er;  // IWYU pragma: keep
}  // namespace PUP
/// \endcond

namespace grmhd::Solutions {

/*!
 * \brief A one-dimensional shock solution for an ideal fluid in Minkowski
 * spacetime
 *
 * This solution consists of a left state for \f$x<0\f$ and a right state for
 * \f$x\ge 0\f$, each with constant fluid variables. The interface between these
 * states moves with the shock speed \f$\mu\f$ as described in
 * \cite Komissarov1999.
 *
 * \note We do not currently support 1D RMHD, so this class provides a 3D
 * solution with \f$x\f$-dependence only. Therefore the computational domain can
 * be represented by a single element with periodic boundary conditions in the
 * \f$y\f$ and \f$z\f$ directions.
 */
class KomissarovShock
    : public evolution::initial_data::InitialData,
      public AnalyticSolution,
      public hydro::TemperatureInitialization<KomissarovShock>,
      public MarkAsAnalyticSolution {
 public:
  using equation_of_state_type = EquationsOfState::IdealFluid<true>;

  struct AdiabaticIndex {
    using type = double;
    static constexpr Options::String help = {
        "The adiabatic index of the ideal fluid"};
    static type lower_bound() { return 1.0; }
  };
  struct LeftRestMassDensity {
    using type = double;
    static std::string name() { return "LeftDensity"; };
    static constexpr Options::String help = {
        "Fluid rest mass density in the left half-domain"};
    static type lower_bound() { return 0.0; }
  };
  struct RightRestMassDensity {
    using type = double;
    static std::string name() { return "RightDensity"; };
    static constexpr Options::String help = {
        "Fluid rest mass density in the right half-domain"};
    static type lower_bound() { return 0.0; }
  };

  struct LeftElectronFraction {
    using type = double;
    static std::string name() { return "LeftElectronFraction"; };
    static constexpr Options::String help = {
        "Fluid electron fraction in the left half-domain"};
    static type lower_bound() { return 0.0; }
    static type upper_bound() { return 1.0; }
  };
  struct RightElectronFraction {
    using type = double;
    static std::string name() { return "RightElectronFraction"; };
    static constexpr Options::String help = {
        "Fluid electron fraction in the right half-domain"};
    static type lower_bound() { return 0.0; }
    static type upper_bound() { return 1.0; }
  };
  struct LeftPressure {
    using type = double;
    static constexpr Options::String help = {
        "Fluid pressure in the left half-domain"};
    static type lower_bound() { return 0.0; }
  };
  struct RightPressure {
    using type = double;
    static constexpr Options::String help = {
        "Fluid pressure in the right half-domain"};
    static type lower_bound() { return 0.0; }
  };
  struct LeftSpatialVelocity {
    using type = std::array<double, 3>;
    static std::string name() { return "LeftVelocity"; };
    static constexpr Options::String help = {
        "Fluid spatial velocity in the left half-domain"};
  };
  struct RightSpatialVelocity {
    using type = std::array<double, 3>;
    static std::string name() { return "RightVelocity"; };
    static constexpr Options::String help = {
        "Fluid spatial velocity in the right half-domain"};
  };
  struct LeftMagneticField {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {
        "Magnetic field in the left half-domain"};
  };
  struct RightMagneticField {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {
        "Magnetic field in the right half-domain"};
  };
  struct ShockSpeed {
    using type = double;
    static constexpr Options::String help = {"Propagation speed of the shock"};
  };

  using options =
      tmpl::list<AdiabaticIndex, LeftRestMassDensity, RightRestMassDensity,
                 LeftElectronFraction, RightElectronFraction, LeftPressure,
                 RightPressure, LeftSpatialVelocity, RightSpatialVelocity,
                 LeftMagneticField, RightMagneticField, ShockSpeed>;

  static constexpr Options::String help = {
      "Analytic initial data for a Komissarov shock test. The fluid variables "
      "are set homogeneously on either half of the domain left and right of "
      "x=0."};

  KomissarovShock() = default;
  KomissarovShock(const KomissarovShock& /*rhs*/) = default;
  KomissarovShock& operator=(const KomissarovShock& /*rhs*/) = default;
  KomissarovShock(KomissarovShock&& /*rhs*/) = default;
  KomissarovShock& operator=(KomissarovShock&& /*rhs*/) = default;
  ~KomissarovShock() override = default;

  KomissarovShock(double adiabatic_index, double left_rest_mass_density,
                  double right_rest_mass_density, double left_electron_fraction,
                  double right_electron_fraction, double left_pressure,
                  double right_pressure,
                  const std::array<double, 3>& left_spatial_velocity,
                  const std::array<double, 3>& right_spatial_velocity,
                  const std::array<double, 3>& left_magnetic_field,
                  const std::array<double, 3>& right_magnetic_field,
                  double shock_speed);

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit KomissarovShock(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(KomissarovShock);
  /// \endcond

  /// @{
  /// Retrieve the GRMHD variables at a given position.
  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<hydro::Tags::RestMassDensity<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::RestMassDensity<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<hydro::Tags::ElectronFraction<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::ElectronFraction<DataType>>;

  template <typename DataType>
  auto variables(
      const tnsr::I<DataType, 3>& x, double t,
      tmpl::list<hydro::Tags::SpecificInternalEnergy<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::SpecificInternalEnergy<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<hydro::Tags::Pressure<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::Pressure<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<hydro::Tags::SpatialVelocity<DataType, 3>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::SpatialVelocity<DataType, 3>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<hydro::Tags::MagneticField<DataType, 3>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::MagneticField<DataType, 3>>;

  template <typename DataType>
  auto variables(
      const tnsr::I<DataType, 3>& x, double t,
      tmpl::list<hydro::Tags::DivergenceCleaningField<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::DivergenceCleaningField<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<hydro::Tags::LorentzFactor<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::LorentzFactor<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<hydro::Tags::SpecificEnthalpy<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::SpecificEnthalpy<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<hydro::Tags::Temperature<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::Temperature<DataType>> {
    return TemperatureInitialization::variables(
        x, t, tmpl::list<hydro::Tags::Temperature<DataType>>{});
  }
  /// @}

  /// Retrieve a collection of hydro variables at `(x, t)`
  template <typename DataType, typename Tag1, typename Tag2, typename... Tags>
  tuples::TaggedTuple<Tag1, Tag2, Tags...> variables(
      const tnsr::I<DataType, 3>& x, double t,
      tmpl::list<Tag1, Tag2, Tags...> /*meta*/) const {
    return {tuples::get<Tag1>(variables(x, t, tmpl::list<Tag1>{})),
            tuples::get<Tag2>(variables(x, t, tmpl::list<Tag2>{})),
            tuples::get<Tags>(variables(x, t, tmpl::list<Tags>{}))...};
  }

  /// Retrieve the metric variables
  template <typename DataType, typename Tag,
            Requires<tmpl::list_contains_v<
                gr::analytic_solution_tags<3, DataType>, Tag>> = nullptr>
  tuples::TaggedTuple<Tag> variables(const tnsr::I<DataType, 3>& x, double t,
                                     tmpl::list<Tag> /*meta*/) const {
    return background_spacetime_.variables(x, t, tmpl::list<Tag>{});
  }

  const EquationsOfState::IdealFluid<true>& equation_of_state() const {
    return equation_of_state_;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) override;

 protected:
  EquationsOfState::IdealFluid<true> equation_of_state_{};
  gr::Solutions::Minkowski<3> background_spacetime_{};

  double adiabatic_index_ = std::numeric_limits<double>::signaling_NaN();
  double left_rest_mass_density_ = std::numeric_limits<double>::signaling_NaN();
  double right_rest_mass_density_ =
      std::numeric_limits<double>::signaling_NaN();
  double left_electron_fraction_ = std::numeric_limits<double>::signaling_NaN();
  double right_electron_fraction_ =
      std::numeric_limits<double>::signaling_NaN();
  double left_pressure_ = std::numeric_limits<double>::signaling_NaN();
  double right_pressure_ = std::numeric_limits<double>::signaling_NaN();
  std::array<double, 3> left_spatial_velocity_{
      {std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN()}};
  std::array<double, 3> right_spatial_velocity_{
      {std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN()}};
  std::array<double, 3> left_magnetic_field_{
      {std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN()}};
  std::array<double, 3> right_magnetic_field_{
      {std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN()}};
  double shock_speed_ = std::numeric_limits<double>::signaling_NaN();

  friend bool operator==(const KomissarovShock& lhs,
                         const KomissarovShock& rhs);

  friend bool operator!=(const KomissarovShock& lhs,
                         const KomissarovShock& rhs);
};

}  // namespace grmhd::Solutions
