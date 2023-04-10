
// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wredundant-decls"
#include "/usr/local/google-benchmark/1.2/include/benchmark/benchmark.h"
#pragma GCC diagnostic pop
#include <charm++.h>
#include <string>
#include <vector>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/CoordinateMaps/Affine.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/ProductMaps.hpp"
#include "Domain/CoordinateMaps/ProductMaps.tpp"
#include "Domain/Structure/Element.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.tpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/Factory.hpp"
#include "PointwiseFunctions/MathFunctions/PowX.hpp"

// Charm looks for this function but since we build without a main function or
// main module we just have it be empty
extern "C" void CkRegisterMainModule(void) {}

// This file is an example of how to do microbenchmark with Google Benchmark
// https://github.com/google/benchmark
// For two examples in different anonymous namespaces

namespace {
// Benchmark of push_back() in std::vector, following Chandler Carruth's talk
// at CppCon in 2015,
// https://www.youtube.com/watch?v=nXaxk27zwlk

// void bench_create(benchmark::State &state) {
//  while (state.KeepRunning()) {
//    std::vector<int> v;
//    benchmark::DoNotOptimize(&v);
//    static_cast<void>(v);
//  }
// }
// BENCHMARK(bench_create);

// void bench_reserve(benchmark::State &state) {
//  while (state.KeepRunning()) {
//    std::vector<int> v;
//    v.reserve(1);
//    benchmark::DoNotOptimize(v.data());
//  }
// }
// BENCHMARK(bench_reserve);

// void bench_push_back(benchmark::State &state) {
//  while (state.KeepRunning()) {
//    std::vector<int> v;
//    v.reserve(1);
//    benchmark::DoNotOptimize(v.data());
//    v.push_back(42);
//    benchmark::ClobberMemory();
//  }
// }
// BENCHMARK(bench_push_back);
}  // namespace

namespace {
// In this anonymous namespace is an example of microbenchmarking the
// all_gradient routine for the GH system

template <size_t Dim>
struct Kappa : db::SimpleTag {
  using type = tnsr::abb<DataVector, Dim, Frame::Grid>;
};
template <size_t Dim>
struct Psi : db::SimpleTag {
  using type = tnsr::aa<DataVector, Dim, Frame::Grid>;
};

// clang-tidy: don't pass be non-const reference
void bench_all_gradient(benchmark::State& state) {  // NOLINT
  constexpr const size_t pts_1d = 4;
  constexpr const size_t Dim = 3;
  const Mesh<Dim> mesh{pts_1d, Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
  domain::CoordinateMaps::Affine map1d(-1.0, 1.0, -1.0, 1.0);
  using Map3d =
      domain::CoordinateMaps::ProductOf3Maps<domain::CoordinateMaps::Affine,
                                             domain::CoordinateMaps::Affine,
                                             domain::CoordinateMaps::Affine>;
  domain::CoordinateMap<Frame::ElementLogical, Frame::Grid, Map3d> map(
      Map3d{map1d, map1d, map1d});

  using VarTags = tmpl::list<Kappa<Dim>, Psi<Dim>>;
  const InverseJacobian<DataVector, Dim, Frame::ElementLogical, Frame::Grid>
      inv_jac = map.inv_jacobian(logical_coordinates(mesh));
  const auto grid_coords = map(logical_coordinates(mesh));
  Variables<VarTags> vars(mesh.number_of_grid_points(), 0.0);

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(partial_derivatives<VarTags>(vars, mesh, inv_jac));
  }
}
// BENCHMARK(bench_all_gradient);  // NOLINT

// BENCHMARK(bench_all_gradient);  // NOLINT

void bench_polytrope_pressure_from_density(benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const auto eos = EquationsOfState::PolytropicFluid<true>{100.0, 2.1};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(eos.pressure_from_density(density));
  }
}
BENCHMARK(bench_polytrope_pressure_from_density);

void bench_polytrope_internal_energy_from_density(
    benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const auto eos = EquationsOfState::PolytropicFluid<true>{100.0, 2.1};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(
        eos.specific_internal_energy_from_density(density));
  }
}
BENCHMARK(bench_polytrope_internal_energy_from_density);

void bench_spectral_pressure_from_density(benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const auto eos = EquationsOfState::Spectral{
      8.2235e-5, 2.5632e-7, {1.35692, 0.0, 0.9297, -0.2523}, 1.5e-3};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(eos.pressure_from_density(density));
  }
}
BENCHMARK(bench_spectral_pressure_from_density);
void bench_spectral_internal_energy_from_density(
    benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const auto eos = EquationsOfState::Spectral{
      8.2235e-5, 2.5632e-7, {1.35692, 0.0, 0.9297, -0.2523}, 1.5e-3};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(
        eos.specific_internal_energy_from_density(density));
  }
}
BENCHMARK(bench_spectral_internal_energy_from_density);

void bench_enthalpy_pressure_from_density(benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const double reference_density{0.0004533804669935759};
  const double max_density{0.0031736632689550316};
  const double min_density{0.0004533804669935759};
  const double trig_scale{0.0};

  const std::vector<double> polynomial_coefficients{
      1.0906760933987152,     0.09067609339871513,   0.045338046699357565,
      0.015112682233119188,   0.003778170558279797,  0.0007556341116559594,
      0.00012593901860932658, 1.7991288372760937e-05};
  const std::vector<double> sin_coefficients{};
  const std::vector<double> cos_coefficients{};
  const EquationsOfState::PolytropicFluid<true> low_density_eos{100.0, 2.0};
  const double transition_delta_eps{0.0};
  const auto eos =
      EquationsOfState::Enthalpy<EquationsOfState::PolytropicFluid<true>>{
          reference_density,
          max_density,
          min_density,
          trig_scale,
          polynomial_coefficients,
          sin_coefficients,
          cos_coefficients,
          low_density_eos,
          transition_delta_eps};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(eos.pressure_from_density(density));
  }
}
BENCHMARK(bench_enthalpy_pressure_from_density);
void bench_enthalpy_internal_energy_from_density(
    benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const double reference_density{0.0004533804669935759};
  const double max_density{0.0031736632689550316};
  const double min_density{0.0004533804669935759};
  const double trig_scale{0.0};

  const std::vector<double> polynomial_coefficients{
      1.0906760933987152,     0.09067609339871513,   0.045338046699357565,
      0.015112682233119188,   0.003778170558279797,  0.0007556341116559594,
      0.00012593901860932658, 1.7991288372760937e-05};
  const std::vector<double> sin_coefficients{};
  const std::vector<double> cos_coefficients{};
  const EquationsOfState::PolytropicFluid<true> low_density_eos{100.0, 2.0};
  const double transition_delta_eps{0.0};
  const auto eos =
      EquationsOfState::Enthalpy<EquationsOfState::PolytropicFluid<true>>{
          reference_density,
          max_density,
          min_density,
          trig_scale,
          polynomial_coefficients,
          sin_coefficients,
          cos_coefficients,
          low_density_eos,
          transition_delta_eps};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(
        eos.specific_internal_energy_from_density(density));
  }
}
BENCHMARK(bench_enthalpy_internal_energy_from_density);

void bench_enthalpy_chi_from_density(benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const double reference_density{0.0004533804669935759};
  const double max_density{0.0031736632689550316};
  const double min_density{0.0004533804669935759};
  const double trig_scale{0.0};

  const std::vector<double> polynomial_coefficients{
      1.0906760933987152,     0.09067609339871513,   0.045338046699357565,
      0.015112682233119188,   0.003778170558279797,  0.0007556341116559594,
      0.00012593901860932658, 1.7991288372760937e-05};
  const std::vector<double> sin_coefficients{};
  const std::vector<double> cos_coefficients{};
  const EquationsOfState::PolytropicFluid<true> low_density_eos{100.0, 2.0};
  const double transition_delta_eps{0.0};
  const auto eos =
      EquationsOfState::Enthalpy<EquationsOfState::PolytropicFluid<true>>{
          reference_density,
          max_density,
          min_density,
          trig_scale,
          polynomial_coefficients,
          sin_coefficients,
          cos_coefficients,
          low_density_eos,
          transition_delta_eps};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(eos.chi_from_density(density));
  }
}
BENCHMARK(bench_enthalpy_chi_from_density);
void bench_enthalpy_sly_pressure_from_density(
    benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const double reference_density{0.00022669023349678794};
  const double max_density{0.0031736632689550316};
  const double min_density{0.0004533804669935759};
  const double trig_scale{1.8428486108259632};

  const std::vector<double> polynomial_coefficients{
      1.0,
      0.02512680820328788,
      7.511888526921824e-26,
      0.006832118818976958,
      0.01663614586347631,
      7.4103163875837e-36,
      4.206931653685116e-49,
      2.6484288063816473e-32,
      1.7304939105067978e-32,
  };
  const std::vector<double> sin_coefficients{0.0012526361570758115,
                                             0.004707127928103492};
  const std::vector<double> cos_coefficients{-0.000578031930149701,
                                             -0.008561072026836546};
  const EquationsOfState::Spectral low_density_eos{
      9.067609339871518e-05,
      3.2387276859504656e-07,
      {1.35, 0.0, -0.7727981253248064, 1.1544540103223992},
      0.0004533804669935759};
  const double transition_delta_eps{0.0};
  const auto eos = EquationsOfState::Enthalpy<EquationsOfState::Spectral>{
      reference_density,
      max_density,
      min_density,
      trig_scale,
      polynomial_coefficients,
      sin_coefficients,
      cos_coefficients,
      low_density_eos,
      transition_delta_eps};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(eos.pressure_from_density(density));
  }
}
BENCHMARK(bench_enthalpy_sly_pressure_from_density);
void bench_enthalpy_sly_internal_energy_from_density(
    benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const double reference_density{0.00022669023349678794};
  const double max_density{0.0031736632689550316};
  const double min_density{0.0004533804669935759};
  const double trig_scale{1.8428486108259632};

  const std::vector<double> polynomial_coefficients{1.0,
                                                    0.02512680820328788,
                                                    7.511888526921824e-26,
                                                    0.006832118818976958,
                                                    0.01663614586347631,
                                                    7.4103163875837e-36,
                                                    4.206931653685116e-49,
                                                    2.6484288063816473e-32,
                                                    1.7304939105067978e-32};
  const std::vector<double> sin_coefficients{0.0012526361570758115,
                                             0.004707127928103492};
  const std::vector<double> cos_coefficients{-0.000578031930149701,
                                             -0.008561072026836546};
  const EquationsOfState::Spectral low_density_eos{
      9.067609339871518e-05,
      3.2387276859504656e-07,
      {1.35, 0.0, -0.7727981253248064, 1.1544540103223992},
      0.0004533804669935759};
  const double transition_delta_eps{0.0};
  const auto eos = EquationsOfState::Enthalpy<EquationsOfState::Spectral>{
      reference_density,
      max_density,
      min_density,
      trig_scale,
      polynomial_coefficients,
      sin_coefficients,
      cos_coefficients,
      low_density_eos,
      transition_delta_eps};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(
        eos.specific_internal_energy_from_density(density));
  }
}
BENCHMARK(bench_enthalpy_sly_internal_energy_from_density);

// This is not actually dbhf, just has a comparable numebr of parameters
void bench_enthalpy_sly_more_pressure_from_density(
    benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const double reference_density{0.00022669023349678794};
  const double max_density{0.0031736632689550316};
  const double min_density{0.0004533804669935759};

  const double trig_scale{0.9759906329155857};

  const std::vector<double> polynomial_coefficients{1.0,
                                                    0.02512680820328788,
                                                    7.511888526921824e-26,
                                                    0.006832118818976958,
                                                    0.01663614586347631,
                                                    7.4103163875837e-36,
                                                    4.206931653685116e-49,
                                                    2.6484288063816473e-32,
                                                    1.7304939105067978e-32};
  const std::vector<double> sin_coefficients{
      -0.2910880174221773, -0.8043646713945052, 0.5415898535044361,
      0.014203340496991165, -0.02097610070937428

  };
  const std::vector<double> cos_coefficients{
      0.8413309497450576, -0.6624885514597745, -0.28831885237845367,
      0.18072953664941335, -0.006075287459066015};
  const EquationsOfState::Spectral low_density_eos{
      9.067609339871518e-05,
      2.9420389962289594e-07,
      {1.35692, 0.0, 1.21113713100009, -0.4058745603468004},
      0.0004533804669935759};
  const double transition_delta_eps{0.0};
  const auto eos = EquationsOfState::Enthalpy<EquationsOfState::Spectral>{
      reference_density,
      max_density,
      min_density,
      trig_scale,
      polynomial_coefficients,
      sin_coefficients,
      cos_coefficients,
      low_density_eos,
      transition_delta_eps};
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(eos.pressure_from_density(density));
  }
}
BENCHMARK(bench_enthalpy_sly_more_pressure_from_density);
void bench_enthalpy_sly_more_internal_energy_from_density(
    benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const double reference_density{0.00022669023349678794};
  const double max_density{0.0031736632689550316};
  const double min_density{0.0004533804669935759};

  const double trig_scale{0.9759906329155857};

  const std::vector<double> polynomial_coefficients{1.0,
                                                    0.02512680820328788,
                                                    7.511888526921824e-26,
                                                    0.006832118818976958,
                                                    0.01663614586347631,
                                                    7.4103163875837e-36,
                                                    4.206931653685116e-49,
                                                    2.6484288063816473e-32,
                                                    1.7304939105067978e-32};
  const std::vector<double> sin_coefficients{
      -0.2910880174221773, -0.8043646713945052, 0.5415898535044361,
      0.014203340496991165, -0.02097610070937428

  };
  const std::vector<double> cos_coefficients{
      0.8413309497450576, -0.6624885514597745, -0.28831885237845367,
      0.18072953664941335, -0.006075287459066015};
  const EquationsOfState::Spectral low_density_eos{
      9.067609339871518e-05,
      2.9420389962289594e-07,
      {1.35692, 0.0, 1.21113713100009, -0.4058745603468004},
      0.0004533804669935759};
  const double transition_delta_eps{0.0};
  const auto eos = EquationsOfState::Enthalpy<EquationsOfState::Spectral>{
      reference_density,
      max_density,
      min_density,
      trig_scale,
      polynomial_coefficients,
      sin_coefficients,
      cos_coefficients,
      low_density_eos,
      transition_delta_eps};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(
        eos.specific_internal_energy_from_density(density));
  }
}
BENCHMARK(bench_enthalpy_sly_more_internal_energy_from_density);
void bench_enthalpy_fast_pressure_from_density(
    benchmark::State& state) {  // NOLINT
  const Scalar<double> density(5e-4);
  const double reference_density{0.0004533804669935759};
  const double max_density{0.0031736632689550316};
  const double min_density{0.0004533804669935759};
  const double trig_scale{0.0};

  const std::vector<double> polynomial_coefficients{1.0, 1.0e-3};
  const std::vector<double> sin_coefficients{};
  const std::vector<double> cos_coefficients{};
  const EquationsOfState::PolytropicFluid<true> low_density_eos{100.0, 2.0};
  const double transition_delta_eps{0.0};
  const auto eos =
      EquationsOfState::Enthalpy<EquationsOfState::PolytropicFluid<true>>{
          reference_density,
          max_density,
          min_density,
          trig_scale,
          polynomial_coefficients,
          sin_coefficients,
          cos_coefficients,
          low_density_eos,
          transition_delta_eps};

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(eos.pressure_from_density(density));
  }
}
BENCHMARK(bench_enthalpy_fast_pressure_from_density);
}  // namespace

// Ignore the warning about an extra ';' because some versions of benchmark
// require it
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
BENCHMARK_MAIN();
#pragma GCC diagnostic pop
