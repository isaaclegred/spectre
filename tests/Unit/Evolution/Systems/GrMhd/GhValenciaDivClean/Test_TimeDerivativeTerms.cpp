// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <memory>
#include <random>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/DampedHarmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/System.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/TimeDerivativeTerms.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/System.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/Tags.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/TimeDerivativeTerms.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/DerivSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfLapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfShift.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

namespace {
template <typename ComputeVolumeTimeDerivativeTerms, size_t Dim,
          typename EvolvedTagList, typename FluxTagList, typename TempTagList,
          typename GradientTagList, typename ArgTagList>
struct ComputeVolumeTimeDerivativeTermsHelper;

template <typename ComputeVolumeTimeDerivativeTerms, size_t Dim,
          typename... EvolvedTags, typename... FluxTags, typename... TempTags,
          typename... GradientTags, typename... ArgTags>
struct ComputeVolumeTimeDerivativeTermsHelper<
    ComputeVolumeTimeDerivativeTerms, Dim, tmpl::list<EvolvedTags...>,
    tmpl::list<FluxTags...>, tmpl::list<TempTags...>,
    tmpl::list<GradientTags...>, tmpl::list<ArgTags...>> {
  template <typename EvolvedVariables, typename FluxVariables,
            typename TemporaryVariables, typename GradientVariables,
            typename ArgumentVariables>
  static void apply(
      const gsl::not_null<EvolvedVariables*> dt_vars_ptr,
      [[maybe_unused]] const gsl::not_null<FluxVariables*> volume_fluxes,
      const gsl::not_null<TemporaryVariables*> temporaries,
      const GradientVariables& partial_derivs,
      const ArgumentVariables& time_derivative_args) {
    ComputeVolumeTimeDerivativeTerms::apply(
        make_not_null(&get<::Tags::dt<EvolvedTags>>(*dt_vars_ptr))...,
        make_not_null(&get<FluxTags>(*volume_fluxes))...,
        make_not_null(&get<TempTags>(*temporaries))...,
        get<GradientTags>(partial_derivs)..., [](const auto& t) -> const auto& {
          if constexpr (tt::is_a_v<std::unique_ptr,
                                   std::decay_t<decltype(t)>>) {
            return *t;
          } else {
            return t;
          }
        }(tuples::get<ArgTags>(time_derivative_args))...);
  }

  template <typename EvolvedVariables, typename FluxVariables,
            typename TemporaryVariables, typename GradientVariables,
            typename ArgumentVariables>
  static void apply_packed(
      const gsl::not_null<EvolvedVariables*> dt_vars_ptr,
      [[maybe_unused]] const gsl::not_null<FluxVariables*> volume_fluxes,
      const gsl::not_null<TemporaryVariables*> temporaries,
      const GradientVariables& partial_derivs,
      const ArgumentVariables& time_derivative_args) {
    ComputeVolumeTimeDerivativeTerms::apply(
        dt_vars_ptr, volume_fluxes, temporaries,
        get<GradientTags>(partial_derivs)..., [](const auto& t) -> const auto& {
          if constexpr (tt::is_a_v<std::unique_ptr,
                                   std::decay_t<decltype(t)>>) {
            return *t;
          } else {
            return t;
          }
        }(tuples::get<ArgTags>(time_derivative_args))...);
  }
};
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GhValenciaDivClean.TimeDerivativeTerms",
    "[Unit][Evolution]") {
  using gh_variables_tags = typename gh::System<3_st>::variables_tag::tags_list;
  using valencia_variables_tags =
      typename grmhd::ValenciaDivClean::System::variables_tag::tags_list;

  using gh_dt_variables_tags =
      grmhd::GhValenciaDivClean::TimeDerivativeTerms::gh_dt_tags;
  using valencia_dt_variables_tags =
      grmhd::GhValenciaDivClean::TimeDerivativeTerms::valencia_dt_tags;
  using dt_variables_type =
      Variables<tmpl::append<gh_dt_variables_tags, valencia_dt_variables_tags>>;

  using gh_flux_tags = tmpl::list<>;
  using valencia_flux_tags =
      grmhd::GhValenciaDivClean::TimeDerivativeTerms::valencia_flux_tags;
  using flux_variables_type =
      Variables<tmpl::append<gh_flux_tags, valencia_flux_tags>>;

  using gh_temp_tags =
      grmhd::GhValenciaDivClean::TimeDerivativeTerms::gh_temp_tags;
  using valencia_temp_tags =
      grmhd::GhValenciaDivClean::TimeDerivativeTerms::valencia_temp_tags;
  using temp_variables_type = Variables<
      typename grmhd::GhValenciaDivClean::TimeDerivativeTerms::temporary_tags>;

  using gh_gradient_tags = tmpl::transform<
      grmhd::GhValenciaDivClean::TimeDerivativeTerms::gh_gradient_tags,
      tmpl::bind<::Tags::deriv, tmpl::_1, tmpl::pin<tmpl::size_t<3_st>>,
                 tmpl::pin<Frame::Inertial>>>;
  using valencia_gradient_tags = tmpl::list<>;
  using gradient_variables_type = Variables<gh_gradient_tags>;

  using gh_arg_tags =
      grmhd::GhValenciaDivClean::TimeDerivativeTerms::gh_arg_tags;
  using valencia_arg_tags =
      grmhd::GhValenciaDivClean::TimeDerivativeTerms::valencia_arg_tags;
  using all_valencia_arg_tags =
      typename grmhd::ValenciaDivClean::TimeDerivativeTerms::argument_tags;
  using SpatialMetricTag = gr::Tags::SpatialMetric<DataVector, 3>;
  using d_SpatialMetricTag =
      ::Tags::deriv<SpatialMetricTag, tmpl::size_t<3>, Frame::Inertial>;
  using arg_variables_type = tuples::tagged_tuple_from_typelist<
      tmpl::remove<tmpl::append<gh_arg_tags, valencia_arg_tags>,
                   gr::Tags::SpatialMetric<DataVector, 3>>>;

  const size_t element_size = 10_st;
  MAKE_GENERATOR(gen);
  std::uniform_real_distribution<> dist(0.1, 1.0);

  dt_variables_type expected_dt_variables{element_size};
  dt_variables_type dt_variables{element_size};
  // Because the tilde_d evolution equation has no source terms, the dt_tilde_d
  // isn't set by the valencia div clean system, so we just set it arbitrarily
  fill_with_random_values(
      make_not_null(&get<::Tags::dt<grmhd::ValenciaDivClean::Tags::TildeD>>(
          expected_dt_variables)),
      make_not_null(&gen), make_not_null(&dist));
  get<::Tags::dt<grmhd::ValenciaDivClean::Tags::TildeD>>(dt_variables) =
      get<::Tags::dt<grmhd::ValenciaDivClean::Tags::TildeD>>(
          expected_dt_variables);

  fill_with_random_values(
      make_not_null(&get<::Tags::dt<grmhd::ValenciaDivClean::Tags::TildeYe>>(
          expected_dt_variables)),
      make_not_null(&gen), make_not_null(&dist));
  get<::Tags::dt<grmhd::ValenciaDivClean::Tags::TildeYe>>(dt_variables) =
      get<::Tags::dt<grmhd::ValenciaDivClean::Tags::TildeYe>>(
          expected_dt_variables);

  flux_variables_type expected_flux_variables{element_size};
  flux_variables_type flux_variables{element_size};

  temp_variables_type temp_variables{element_size};
  temp_variables_type expected_temp_variables{element_size};

  const auto gradient_variables =
      make_with_random_values<gradient_variables_type>(
          make_not_null(&gen), make_not_null(&dist), element_size);
  arg_variables_type arg_variables;
  tmpl::for_each<
      tmpl::remove<tmpl::remove<tmpl::append<gh_arg_tags, valencia_arg_tags>,
                                SpatialMetricTag>,
                   d_SpatialMetricTag>>([&gen, &dist,
                                         &arg_variables](auto tag_v) {
    using tag = typename decltype(tag_v)::type;
    if constexpr (std::is_same_v<
                      typename tag::type,
                      std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>>) {
      tuples::get<tag>(arg_variables) = make_with_random_values<
          typename tnsr::I<DataVector, 3, Frame::Inertial>>(
          make_not_null(&gen), make_not_null(&dist), DataVector{element_size});
    } else if constexpr (tt::is_a_v<Tensor, typename tag::type>) {
      tuples::get<tag>(arg_variables) =
          make_with_random_values<typename tag::type>(make_not_null(&gen),
                                                      make_not_null(&dist),
                                                      DataVector{element_size});
    }
  });
  get<gh::gauges::Tags::GaugeCondition>(arg_variables) =
      std::make_unique<gh::gauges::DampedHarmonic>(
          100., std::array{1.2, 1.5, 1.7}, std::array{2, 4, 6});

  // ensure that the signature of the metric is correct
  {
    auto& metric =
        tuples::get<gr::Tags::SpacetimeMetric<DataVector, 3_st>>(arg_variables);
    get<0, 0>(metric) += -2.0;
    for (size_t i = 0; i < 3; ++i) {
      metric.get(i + 1, i + 1) += 4.0;
      metric.get(i + 1, 0) *= 0.01;
    }
  }

  ComputeVolumeTimeDerivativeTermsHelper<
      gh::TimeDerivative<3_st>, 3_st, gh_variables_tags, gh_flux_tags,
      gh_temp_tags, gh_gradient_tags,
      tmpl::remove<gh_arg_tags, gr::Tags::SpatialMetric<DataVector, 3>>>::
      apply(make_not_null(&expected_dt_variables),
            make_not_null(&expected_flux_variables),
            make_not_null(&expected_temp_variables), gradient_variables,
            arg_variables);
  // appropriately mimic the behavior of the `TimeDerivativeTerms`
  // implementation wherein it uses the temporary tags from the Generalized
  // Harmonic system that are applicable to the spacetime parts of the Valencia
  // system.
  tuples::tagged_tuple_from_typelist<all_valencia_arg_tags>
      all_valencia_argument_variables{};

  tmpl::for_each<
      tmpl::remove<all_valencia_arg_tags,
                   SpatialMetricTag>>([&arg_variables, &expected_temp_variables,
                                       &all_valencia_argument_variables](
                                          const auto tag_v) {
    using tag = typename decltype(tag_v)::type;
    if constexpr (tmpl::list_contains_v<gh_temp_tags, tag>) {
      tuples::get<tag>(all_valencia_argument_variables) =
          get<tag>(expected_temp_variables);
    } else if constexpr (std::is_same_v<
                             tag,
                             ::Tags::deriv<gr::Tags::Lapse<DataVector>,
                                           tmpl::size_t<3>, Frame::Inertial>>) {
      tuples::get<tag>(all_valencia_argument_variables) =
          gh::spatial_deriv_of_lapse(
              get<gr::Tags::Lapse<DataVector>>(expected_temp_variables),
              get<gr::Tags::SpacetimeNormalVector<DataVector, 3>>(
                  expected_temp_variables),
              get<gh::Tags::Phi<DataVector, 3>>(arg_variables));
      get<tag>(expected_temp_variables) =
          tuples::get<tag>(all_valencia_argument_variables);
    } else if constexpr (std::is_same_v<
                             tag,
                             ::Tags::deriv<gr::Tags::Shift<DataVector, 3>,
                                           tmpl::size_t<3>, Frame::Inertial>>) {
      tuples::get<tag>(all_valencia_argument_variables) =
          gh::spatial_deriv_of_shift(
              get<gr::Tags::Lapse<DataVector>>(expected_temp_variables),
              get<gr::Tags::InverseSpacetimeMetric<DataVector, 3>>(
                  expected_temp_variables),
              get<gr::Tags::SpacetimeNormalVector<DataVector, 3>>(
                  expected_temp_variables),
              get<gh::Tags::Phi<DataVector, 3>>(arg_variables));
      get<tag>(expected_temp_variables) =
          tuples::get<tag>(all_valencia_argument_variables);
    } else if constexpr (std::is_same_v<tag, gr::Tags::ExtrinsicCurvature<
                                                 DataVector, 3>>) {
      tuples::get<tag>(all_valencia_argument_variables) =
          gh::extrinsic_curvature(
              get<gr::Tags::SpacetimeNormalVector<DataVector, 3>>(
                  expected_temp_variables),
              get<gh::Tags::Pi<DataVector, 3>>(arg_variables),
              get<gh::Tags::Phi<DataVector, 3>>(arg_variables));
      get<tag>(expected_temp_variables) =
          tuples::get<tag>(all_valencia_argument_variables);
    } else {
      tuples::get<tag>(all_valencia_argument_variables) =
          tuples::get<tag>(arg_variables);
    }
  });
  // Set spatial metric and its derivative
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      make_const_view(
          make_not_null(
              &std::as_const(get<gr::Tags::SpatialMetric<DataVector, 3>>(
                                 all_valencia_argument_variables))
                   .get(i, j)),
          get<gr::Tags::SpacetimeMetric<DataVector, 3>>(arg_variables)
              .get(i + 1, j + 1),
          0, element_size);
    }
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = j; k < 3; ++k) {
        make_const_view(
            make_not_null(&std::as_const(get<d_SpatialMetricTag>(
                                             all_valencia_argument_variables))
                               .get(i, j, k)),
            get<gh::Tags::Phi<DataVector, 3>>(arg_variables)
                .get(i, j + 1, k + 1),
            0, element_size);
      }
    }
  }

  ComputeVolumeTimeDerivativeTermsHelper<
      grmhd::ValenciaDivClean::TimeDerivativeTerms, 3_st,
      valencia_variables_tags, valencia_flux_tags, valencia_temp_tags,
      valencia_gradient_tags,
      all_valencia_arg_tags>::apply(make_not_null(&expected_dt_variables),
                                    make_not_null(&expected_flux_variables),
                                    make_not_null(&expected_temp_variables),
                                    gradient_variables,
                                    all_valencia_argument_variables);

  // compute the stress energy for the expected variables
  grmhd::GhValenciaDivClean::trace_reversed_stress_energy(
      make_not_null(
          &get<grmhd::GhValenciaDivClean::Tags::TraceReversedStressEnergy>(
              expected_temp_variables)),
      make_not_null(&get<grmhd::GhValenciaDivClean::Tags::FourVelocityOneForm>(
          expected_temp_variables)),
      make_not_null(
          &get<grmhd::GhValenciaDivClean::Tags::ComovingMagneticFieldOneForm>(
              expected_temp_variables)),
      tuples::get<hydro::Tags::RestMassDensity<DataVector>>(arg_variables),
      tuples::get<hydro::Tags::SpecificEnthalpy<DataVector>>(arg_variables),
      get<hydro::Tags::SpatialVelocityOneForm<DataVector, 3, Frame::Inertial>>(
          expected_temp_variables),
      get<hydro::Tags::MagneticFieldOneForm<DataVector, 3>>(
          expected_temp_variables),
      get<hydro::Tags::MagneticFieldSquared<DataVector>>(
          expected_temp_variables),
      get<hydro::Tags::MagneticFieldDotSpatialVelocity<DataVector>>(
          expected_temp_variables),
      tuples::get<hydro::Tags::LorentzFactor<DataVector>>(arg_variables),
      get<typename grmhd::ValenciaDivClean::TimeDerivativeTerms::
              OneOverLorentzFactorSquared>(expected_temp_variables),
      tuples::get<hydro::Tags::Pressure<DataVector>>(arg_variables),
      tuples::get<gr::Tags::SpacetimeMetric<DataVector, 3_st>>(arg_variables),
      get<gr::Tags::Shift<DataVector, 3_st>>(expected_temp_variables),
      get<gr::Tags::Lapse<DataVector>>(expected_temp_variables));
  // apply the correction to dt pi for the expected variables
  grmhd::GhValenciaDivClean::add_stress_energy_term_to_dt_pi(
      make_not_null(&get<::Tags::dt<gh::Tags::Pi<DataVector, 3_st>>>(
          expected_dt_variables)),
      get<grmhd::GhValenciaDivClean::Tags::TraceReversedStressEnergy>(
          expected_temp_variables),
      get<gr::Tags::Lapse<DataVector>>(expected_temp_variables));

  ComputeVolumeTimeDerivativeTermsHelper<
      grmhd::GhValenciaDivClean::TimeDerivativeTerms, 3_st,
      tmpl::append<gh_variables_tags, valencia_variables_tags>,
      typename flux_variables_type::tags_list,
      typename temp_variables_type::tags_list,
      typename gradient_variables_type::tags_list,
      tmpl::remove<tmpl::remove<typename arg_variables_type::tags_list,
                                SpatialMetricTag>,
                   d_SpatialMetricTag>>::
      apply_packed(make_not_null(&dt_variables), make_not_null(&flux_variables),
                   make_not_null(&temp_variables), gradient_variables,
                   arg_variables);

  CHECK_VARIABLES_APPROX(dt_variables, expected_dt_variables);
  CHECK_VARIABLES_APPROX(flux_variables, expected_flux_variables);
  CHECK_VARIABLES_APPROX(temp_variables, expected_temp_variables);
}
