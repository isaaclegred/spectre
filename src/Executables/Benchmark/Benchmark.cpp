// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wredundant-decls"
#include <benchmark/benchmark.h>
#pragma GCC diagnostic pop
#include <charm++.h>
#include <string>
#include <vector>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Identity.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/CoordinateMaps/Affine.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/ProductMaps.hpp"
#include "Domain/CoordinateMaps/ProductMaps.tpp"
#include "Domain/Structure/Element.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/ConservativeFromPrimitive.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/KastaunEtAl.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/NewmanHamlin.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/PalenzuelaEtAl.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/PrimitiveFromConservative.hpp"
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
void bench_con_2_prim(benchmark::State& state) {
  EquationsOfState::PolytropicFluid<true> eos(100.0, 2.0);
  Scalar<DataVector> rest_mass_density{
      DataVector{1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3}};
  Scalar<DataVector> electron_fraction{
      DataVector{.1, .1, .1, .1, .1, .1, .1, .1}};
  Scalar<DataVector> specific_internal_energy =
      eos.specific_internal_energy_from_density(rest_mass_density);
  tnsr::I<DataVector, 3, Frame::Inertial> spatial_velocity{
      {DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}};
  tnsr::I<DataVector, 3, Frame::Inertial> magnetic_field{
      {DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}};

  Scalar<DataVector> divergence_cleaning_field{
      DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  Scalar<DataVector> lorentz_factor{
      DataVector{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
  Scalar<DataVector> pressure = eos.pressure_from_density(rest_mass_density);
  Scalar<DataVector> specific_enthalpy =
      Scalar<DataVector>(get(pressure) / get(rest_mass_density) +
                         get(specific_internal_energy) + 1.0);
  tnsr::i<DataVector, 3, Frame::Inertial> tilde_s{
      {DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}};
  tnsr::I<DataVector, 3, Frame::Inertial> tilde_b{
      {DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}};
  Scalar<DataVector> tilde_d{
      DataVector{1e-10, 1e-7, 1e-5, 1e-4, 5e-4, 1e-3, 2e-3}};
  Scalar<DataVector> tilde_ye{DataVector{.1, .1, .1, .1, .1, .1, .1, .1}};
  Scalar<DataVector> tilde_tau{
      DataVector{1.1e-10, 1.2e-7, 1.5e-5, 1.7e-4, 6e-4, 2.1e-3, 6e-3}};
  Scalar<DataVector> tilde_phi{DataVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  const Scalar<DataVector> sqrt_det_spatial_metric{
      DataVector{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
  size_t Dim = 3;
  DataVector used_for_size{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  tnsr::ii<DataVector, 3, Frame::Inertial> spatial_metric =
      make_with_value<tnsr::ii<DataVector, 3, Frame::Inertial>>(used_for_size,
                                                                0.0);
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = i; j < Dim; ++j) {
      spatial_metric.get(i, j) =
          make_with_value<DataVector>(used_for_size, (i + 1.) * (j + 1.));
    }
  }
  tnsr::II<DataVector, 3, Frame::Inertial> inv_spatial_metric =
      make_with_value<tnsr::II<DataVector, 3, Frame::Inertial>>(used_for_size,
                                                                0.0);
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = i; j < Dim; ++j) {
      spatial_metric.get(i, j) =
          make_with_value<DataVector>(used_for_size, (i + 1.) * (j + 1.));
    }
  }
  grmhd::ValenciaDivClean::ConservativeFromPrimitive::apply(
      make_not_null(&tilde_d), make_not_null(&tilde_ye),
      make_not_null(&tilde_tau), make_not_null(&tilde_s),
      make_not_null(&tilde_b), make_not_null(&tilde_phi), rest_mass_density,
      electron_fraction, specific_internal_energy, pressure, spatial_velocity,
      lorentz_factor, magnetic_field, sqrt_det_spatial_metric, spatial_metric,
      divergence_cleaning_field);
  using ordered_list_of_primitive_recovery_schemes = tmpl::list<
      grmhd::ValenciaDivClean::PrimitiveRecoverySchemes::KastaunEtAl,
      grmhd::ValenciaDivClean::PrimitiveRecoverySchemes::NewmanHamlin,
      grmhd::ValenciaDivClean::PrimitiveRecoverySchemes::PalenzuelaEtAl>;
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(
        grmhd::ValenciaDivClean::PrimitiveFromConservative<
            ordered_list_of_primitive_recovery_schemes>::
            apply(make_not_null(&rest_mass_density),
                  make_not_null(&electron_fraction),
                  make_not_null(&specific_internal_energy),
                  make_not_null(&spatial_velocity),
                  make_not_null(&magnetic_field),
                  make_not_null(&divergence_cleaning_field),
                  make_not_null(&lorentz_factor), make_not_null(&pressure),
                  make_not_null(&specific_enthalpy), tilde_d, tilde_ye,
                  tilde_tau, tilde_s, tilde_b, tilde_phi, spatial_metric,
                  inv_spatial_metric, sqrt_det_spatial_metric, eos));
  }
}
// In this anonymous namespace is an example of microbenchmarking the
// all_gradient routine for the GH system

// template <size_t Dim>
// struct Kappa : db::SimpleTag {
//   using type = tnsr::abb<DataVector, Dim, Frame::Grid>;
// };
// template <size_t Dim>
// struct Psi : db::SimpleTag {
//   using type = tnsr::aa<DataVector, Dim, Frame::Grid>;
// };

// // clang-tidy: don't pass be non-const reference
// void bench_all_gradient(benchmark::State& state) {  // NOLINT
//   constexpr const size_t pts_1d = 4;
//   constexpr const size_t Dim = 3;
//   const Mesh<Dim> mesh{pts_1d, Spectral::Basis::Legendre,
//                        Spectral::Quadrature::GaussLobatto};
//   domain::CoordinateMaps::Affine map1d(-1.0, 1.0, -1.0, 1.0);
//   using Map3d =
//       domain::CoordinateMaps::ProductOf3Maps<domain::CoordinateMaps::Affine,
//                                              domain::CoordinateMaps::Affine,
//                                              domain::CoordinateMaps::Affine>;
//   domain::CoordinateMap<Frame::ElementLogical, Frame::Grid, Map3d> map(
//       Map3d{map1d, map1d, map1d});

//   using VarTags = tmpl::list<Kappa<Dim>, Psi<Dim>>;
//   const InverseJacobian<DataVector, Dim, Frame::ElementLogical, Frame::Grid>
//       inv_jac = map.inv_jacobian(logical_coordinates(mesh));
//   const auto grid_coords = map(logical_coordinates(mesh));
//   Variables<VarTags> vars(mesh.number_of_grid_points(), 0.0);

//   while (state.KeepRunning()) {
//     benchmark::DoNotOptimize(partial_derivatives<VarTags>(vars, mesh,
//     inv_jac));
//  }

BENCHMARK(bench_con_2_prim);  // NOLINT
}  // namespace

// Ignore the warning about an extra ';' because some versions of benchmark
// require it
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
BENCHMARK_MAIN();
#pragma GCC diagnostic pop
