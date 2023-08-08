// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/DhGaugeParameters.hpp"

#include <array>
#include <pup.h>
#include <pup_stl.h>

gh::gauges::DhGaugeParameters<true>::DhGaugeParameters(
    const double start, const double window, const double width,
    const std::array<double, 3>& amps, const std::array<int, 3>& exps)
    : rollon_start(start),
      rollon_window(window),
      spatial_decay_width(width),
      amplitudes(amps),
      exponents(exps) {}

gh::gauges::DhGaugeParameters<false>::DhGaugeParameters(
    const double width, const std::array<double, 3>& amps,
    const std::array<int, 3>& exps)
    : spatial_decay_width(width), amplitudes(amps), exponents(exps) {}

// NOLINTNEXTLINE(google-runtime-references)
void gh::gauges::DhGaugeParameters<true>::pup(PUP::er& p) {
  p | rollon_start;
  p | rollon_window;
  p | spatial_decay_width;
  p | amplitudes;
  p | exponents;
}

// NOLINTNEXTLINE(google-runtime-references)
void gh::gauges::DhGaugeParameters<false>::pup(PUP::er& p) {
  p | spatial_decay_width;
  p | amplitudes;
  p | exponents;
}
