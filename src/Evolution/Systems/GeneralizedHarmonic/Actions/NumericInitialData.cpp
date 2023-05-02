// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Actions/NumericInitialData.hpp"

#include <boost/functional/hash.hpp>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Phi.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Pi.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/TimeDerivativeOfSpatialMetric.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/PrettyType.hpp"

namespace gh {

void initial_gh_variables_from_adm(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
    const gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const tnsr::ii<DataVector, 3>& extrinsic_curvature, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inv_jacobian) {
  // Assemble spacetime metric from 3+1 vars
  gr::spacetime_metric(spacetime_metric, lapse, shift, spatial_metric);

  // Compute Phi from numerical derivative of the spacetime metric so it
  // satisfies the 3-index constraint
  partial_derivative(phi, *spacetime_metric, mesh, inv_jacobian);

  // Compute Pi by choosing dt_lapse = 0 and dt_shift = 0 (for now).
  // The mutator `SetPiFromGauge` should be combined with this.
  const auto deriv_spatial_metric =
      tenex::evaluate<ti::i, ti::j, ti::k>((*phi)(ti::i, ti::j, ti::k));
  const auto deriv_shift = partial_derivative(shift, mesh, inv_jacobian);
  const auto dt_lapse = make_with_value<Scalar<DataVector>>(lapse, 0.);
  const auto dt_shift = make_with_value<tnsr::I<DataVector, 3>>(shift, 0.);
  const auto dt_spatial_metric = gr::time_derivative_of_spatial_metric(
      lapse, shift, deriv_shift, spatial_metric, deriv_spatial_metric,
      extrinsic_curvature);
  gh::pi(pi, lapse, dt_lapse, shift, dt_shift, spatial_metric,
         dt_spatial_metric, *phi);
}

NumericInitialData::NumericInitialData(
    std::string file_glob, std::string subfile_name,
    std::variant<double, importers::ObservationSelector> observation_value,
    bool enable_interpolation, std::variant<AdmVars, GhVars> selected_variables)
    : importer_options_(std::move(file_glob), std::move(subfile_name),
                        observation_value, enable_interpolation),
      selected_variables_(std::move(selected_variables)) {}

NumericInitialData::NumericInitialData(CkMigrateMessage* msg)
    : InitialData(msg) {}

PUP::able::PUP_ID NumericInitialData::my_PUP_ID = 0;

size_t NumericInitialData::volume_data_id() const {
  size_t hash = 0;
  boost::hash_combine(hash, pretty_type::get_name<NumericInitialData>());
  boost::hash_combine(hash,
                      get<importers::OptionTags::FileGlob>(importer_options_));
  boost::hash_combine(hash,
                      get<importers::OptionTags::Subgroup>(importer_options_));
  return hash;
}

void NumericInitialData::pup(PUP::er& p) {
  p | importer_options_;
  p | selected_variables_;
}

bool operator==(const NumericInitialData& lhs, const NumericInitialData& rhs) {
  return lhs.importer_options_ == rhs.importer_options_ and
         lhs.selected_variables_ == rhs.selected_variables_;
}

}  // namespace gh
