// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/SpinWeighted.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/BouncingBlackHole.hpp"
#include "Evolution/Systems/Cce/BoundaryData.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Helpers/Evolution/Systems/Cce/AnalyticSolutions/AnalyticDataHelpers.hpp"
#include "Helpers/Evolution/Systems/Cce/BoundaryTestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/SwshCollocation.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfLapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfShift.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivOfLapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivOfShift.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivOfSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivativeOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "Utilities/Gsl.hpp"

namespace Cce {
namespace {
tnsr::A<DataVector, 3> bouncing_bh_mapped_coordinates(
    const tnsr::A<DataVector, 3>& unmapped_coordinates, const double amplitude,
    const double period) {
  tnsr::A<DataVector, 3> result{get<0>(unmapped_coordinates).size()};
  get<0>(result) = get<0>(unmapped_coordinates);
  get<1>(result) =
      get<1>(unmapped_coordinates) +
      amplitude *
          pow<4>(sin(2.0 * M_PI * get<0>(unmapped_coordinates) / period));
  get<2>(result) = get<2>(unmapped_coordinates);
  get<3>(result) = get<3>(unmapped_coordinates);
  return result;
}

tnsr::aB<DataVector, 3> bouncing_bh_mapped_jacobian(
    const tnsr::A<DataVector, 3>& unmapped_coordinates, const double amplitude,
    const double period) {
  tnsr::aB<DataVector, 3> result{get<0>(unmapped_coordinates).size(), 0.0};
  get<0, 0>(result) = 1.0;
  get<0, 1>(result) =
      4.0 * amplitude *
      pow<3>(sin(2.0 * M_PI * get<0>(unmapped_coordinates) / period)) *
      cos(2.0 * M_PI * get<0>(unmapped_coordinates) / period) * 2.0 * M_PI /
      period;
  get<1, 1>(result) = 1.0;
  get<2, 2>(result) = 1.0;
  get<3, 3>(result) = 1.0;
  return result;
}

tnsr::aaB<DataVector, 3> bouncing_bh_mapped_hessian(
    const tnsr::A<DataVector, 3>& unmapped_coordinates, const double amplitude,
    const double period) {
  tnsr::aaB<DataVector, 3> result{get<0>(unmapped_coordinates).size(), 0.0};
  get<0, 0, 1>(result) =
      8.0 * amplitude * square(M_PI / period) *
      (cos(4.0 * M_PI * get<0>(unmapped_coordinates) / period) -
       cos(8.0 * M_PI * get<0>(unmapped_coordinates) / period));
  return result;
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.Cce.BouncingBlackHole",
                  "[Unit][Cce]") {
  MAKE_GENERATOR(gen);
  UniformCustomDistribution<double> large_parameter_dist{25.0, 100.0};
  UniformCustomDistribution<double> small_parameter_dist{0.5, 2.0};

  const double amplitude = small_parameter_dist(gen);
  const double period = large_parameter_dist(gen);
  const double mass = small_parameter_dist(gen);

  const double extraction_radius = large_parameter_dist(gen);
  const size_t l_max = 12;
  const double time = large_parameter_dist(gen);

  tnsr::A<DataVector, 3> unmapped_coordinates{
      Spectral::Swsh::number_of_swsh_collocation_points(l_max)};
  const auto& collocation_metadata =
      Spectral::Swsh::cached_collocation_metadata<
          Spectral::Swsh::ComplexRepresentation::Interleaved>(l_max);
  for (const auto collocation_point : collocation_metadata) {
    get<1>(unmapped_coordinates)[collocation_point.offset] =
        extraction_radius * cos(collocation_point.phi) *
        sin(collocation_point.theta);
    get<2>(unmapped_coordinates)[collocation_point.offset] =
        extraction_radius * sin(collocation_point.phi) *
        sin(collocation_point.theta);
    get<3>(unmapped_coordinates)[collocation_point.offset] =
        extraction_radius * cos(collocation_point.theta);
  }
  get<0>(unmapped_coordinates) = time;

  const auto mapped_coordinates =
      bouncing_bh_mapped_coordinates(unmapped_coordinates, amplitude, period);

  const auto mapped_radius = sqrt(square(get<1>(mapped_coordinates)) +
                                  square(get<2>(mapped_coordinates)) +
                                  square(get<3>(mapped_coordinates)));

  tnsr::aa<DataVector, 3> schwarzschild_metric_at_mapped_coordinates{
      Spectral::Swsh::number_of_swsh_collocation_points(l_max)};
  get<0, 0>(schwarzschild_metric_at_mapped_coordinates) =
      2.0 * mass / mapped_radius - 1.0;
  for (size_t i = 0; i < 3; ++i) {
    schwarzschild_metric_at_mapped_coordinates.get(0, i + 1) =
        2.0 * mass / mapped_radius *
        (mapped_coordinates.get(i + 1) / mapped_radius);
    schwarzschild_metric_at_mapped_coordinates.get(i + 1, i + 1) =
        2.0 * mass * square(mapped_coordinates.get(i + 1)) /
            pow<3>(mapped_radius) +
        1.0;
    for (size_t j = 0; j < i; ++j) {
      schwarzschild_metric_at_mapped_coordinates.get(i + 1, j + 1) =
          2.0 * mass * mapped_coordinates.get(i + 1) *
          mapped_coordinates.get(j + 1) / pow<3>(mapped_radius);
    }
  }

  // Apply the Jacobian to determine the down-index metric in the new coordinate
  // system.
  const auto mapped_jacobian =
      bouncing_bh_mapped_jacobian(unmapped_coordinates, amplitude, period);
  tnsr::aa<DataVector, 3> mapped_schwarzschild_metric{
      Spectral::Swsh::number_of_swsh_collocation_points(l_max), 0.0};

  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          mapped_schwarzschild_metric.get(a, b) +=
              mapped_jacobian.get(a, c) * mapped_jacobian.get(b, d) *
              schwarzschild_metric_at_mapped_coordinates.get(c, d);
        }
      }
    }
  }

  tnsr::abb<DataVector, 3> d_schwarzschild_metric_at_mapped_coordinates{
      Spectral::Swsh::number_of_swsh_collocation_points(l_max), 0.0};
  // vanishing time derivative, so we just skip it in the loop
  for (size_t k = 0; k < 3; ++k) {
    d_schwarzschild_metric_at_mapped_coordinates.get(k + 1, 0, 0) =
        -2.0 * mass * mapped_coordinates.get(k + 1) / pow<3>(mapped_radius);
    for (size_t i = 0; i < 3; ++i) {
      if (k == i) {
        d_schwarzschild_metric_at_mapped_coordinates.get(k + 1, 0, i + 1) =
            2.0 * mass / square(mapped_radius) -
            4.0 * mass * square(mapped_coordinates.get(k + 1)) /
                pow<4>(mapped_radius);
      } else {
        d_schwarzschild_metric_at_mapped_coordinates.get(k + 1, 0, i + 1) =
            -4.0 * mass * mapped_coordinates.get(k + 1) *
            mapped_coordinates.get(i + 1) / pow<4>(mapped_radius);
      }
      for (size_t j = i; j < 3; ++j) {
        if (k == i) {
          if (i == j) {
            d_schwarzschild_metric_at_mapped_coordinates.get(k + 1, j + 1,
                                                             i + 1) =
                4.0 * mass * mapped_coordinates.get(j + 1) /
                    pow<3>(mapped_radius) -
                6.0 * mass * mapped_coordinates.get(k + 1) *
                    mapped_coordinates.get(j + 1) *
                    mapped_coordinates.get(i + 1) / pow<5>(mapped_radius);
          } else {
            d_schwarzschild_metric_at_mapped_coordinates.get(k + 1, j + 1,
                                                             i + 1) =
                2.0 * mass * mapped_coordinates.get(j + 1) /
                    pow<3>(mapped_radius) -
                6.0 * mass * mapped_coordinates.get(k + 1) *
                    mapped_coordinates.get(j + 1) *
                    mapped_coordinates.get(i + 1) / pow<5>(mapped_radius);
          }
        } else if (k == j) {
          d_schwarzschild_metric_at_mapped_coordinates.get(k + 1, j + 1,
                                                           i + 1) =
              2.0 * mass * mapped_coordinates.get(i + 1) /
                  pow<3>(mapped_radius) -
              6.0 * mass * mapped_coordinates.get(k + 1) *
                  mapped_coordinates.get(j + 1) *
                  mapped_coordinates.get(i + 1) / pow<5>(mapped_radius);
        } else {
          d_schwarzschild_metric_at_mapped_coordinates.get(k + 1, j + 1,
                                                           i + 1) =
              -6.0 * mass * mapped_coordinates.get(k + 1) *
              mapped_coordinates.get(j + 1) * mapped_coordinates.get(i + 1) /
              pow<5>(mapped_radius);
        }
      }
    }
  }

  const auto mapped_hessian =
      bouncing_bh_mapped_hessian(unmapped_coordinates, amplitude, period);

  tnsr::abb<DataVector, 3> mapped_d_schwarzschild_metric{
      Spectral::Swsh::number_of_swsh_collocation_points(l_max), 0.0};

  for (size_t c = 0; c < 4; ++c) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        for (size_t d = 0; d < 4; ++d) {
          for (size_t e = 0; e < 4; ++e) {
            mapped_d_schwarzschild_metric.get(c, a, b) +=
                mapped_hessian.get(a, c, d) *
                    schwarzschild_metric_at_mapped_coordinates.get(d, e) *
                    mapped_jacobian.get(b, e) +
                mapped_jacobian.get(a, d) * mapped_hessian.get(b, c, e) *
                    schwarzschild_metric_at_mapped_coordinates.get(d, e);
            for (size_t f = 1; f < 4; ++f) {
              mapped_d_schwarzschild_metric.get(c, a, b) +=
                  mapped_jacobian.get(a, d) * mapped_jacobian.get(b, e) *
                  mapped_jacobian.get(c, f) *
                  d_schwarzschild_metric_at_mapped_coordinates.get(f, d, e);
            }
          }
        }
      }
    }
  }
  // Now we test these tensors against the corresponding tensors from the
  // generator
  const auto analytic_solution =
      Solutions::BouncingBlackHole{amplitude, extraction_radius, mass, period};
  const auto boundary_tuple = analytic_solution.variables(
      l_max, time, Solutions::BouncingBlackHole::tags{});
  const auto& spacetime_metric =
      get<gr::Tags::SpacetimeMetric<DataVector, 3>>(boundary_tuple);
  const auto& dt_spacetime_metric =
      get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, 3>>>(boundary_tuple);
  const auto& d_spacetime_metric =
      get<gh::Tags::Phi<DataVector, 3>>(boundary_tuple);

  const auto serialized_and_deserialized_analytic_solution =
      serialize_and_deserialize(analytic_solution);
  const auto boundary_tuple_from_serialized =
      serialized_and_deserialized_analytic_solution.variables(
          l_max, time, tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3>>{});
  const auto& spacetime_metric_from_serialized =
      get<gr::Tags::SpacetimeMetric<DataVector, 3>>(
          boundary_tuple_from_serialized);
  CHECK_ITERABLE_APPROX(spacetime_metric_from_serialized, spacetime_metric);

  CHECK_ITERABLE_APPROX(spacetime_metric, mapped_schwarzschild_metric);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      CHECK_ITERABLE_APPROX(dt_spacetime_metric.get(a, b),
                            mapped_d_schwarzschild_metric.get(0, a, b));
    }
  }

  for (size_t i = 0; i < 3; ++i) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        CHECK_ITERABLE_APPROX(d_spacetime_metric.get(i, a, b),
                              mapped_d_schwarzschild_metric.get(i + 1, a, b));
      }
    }
  }

  // check the coordinate values from the analytic solution (actually testing
  // functionality from the abstract base class `WorldtubeData`).
  const auto& cartesian_coordinates =
      get<Tags::CauchyCartesianCoords>(boundary_tuple);
  const auto& dr_cartesian_coordinates =
      get<Tags::Dr<Tags::CauchyCartesianCoords>>(boundary_tuple);
  const DataVector radius_squared_from_coordinates =
      square(get<0>(cartesian_coordinates)) +
      square(get<1>(cartesian_coordinates)) +
      square(get<2>(cartesian_coordinates));
  for (const auto& val : radius_squared_from_coordinates) {
    CHECK(approx(val) == square(extraction_radius));
  }
  const DataVector radius_from_coordinates =
      (get<0>(cartesian_coordinates) * get<0>(dr_cartesian_coordinates) +
       get<1>(cartesian_coordinates) * get<1>(dr_cartesian_coordinates) +
       get<2>(cartesian_coordinates) * get<2>(dr_cartesian_coordinates));
  for (const auto& val : radius_from_coordinates) {
    CHECK(approx(val) == extraction_radius);
  }

  Solutions::TestHelpers::check_adm_metric_quantities(
      boundary_tuple, mapped_schwarzschild_metric, dt_spacetime_metric,
      d_spacetime_metric);
  Solutions::TestHelpers::test_initialize_j(
      l_max, 5_st, extraction_radius, time,
      std::make_unique<InitializeJ::InverseCubic<false>>(),
      analytic_solution.get_clone());
}
}  // namespace Cce
