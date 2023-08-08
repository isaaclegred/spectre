// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/Ccz4/ATilde.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace Ccz4 {
template <typename DataType, size_t Dim, typename Frame>
void a_tilde(const gsl::not_null<tnsr::ii<DataType, Dim, Frame>*> result,
             const gsl::not_null<Scalar<DataType>*> buffer,
             const Scalar<DataType>& conformal_factor_squared,
             const tnsr::ii<DataType, Dim, Frame>& spatial_metric,
             const tnsr::ii<DataType, Dim, Frame>& extrinsic_curvature,
             const Scalar<DataType>& trace_extrinsic_curvature) {
  ::tenex::evaluate(buffer, trace_extrinsic_curvature() / 3.0);
  ::tenex::evaluate<ti::i, ti::j>(
      result, conformal_factor_squared() *
                  (extrinsic_curvature(ti::i, ti::j) -
                   (*buffer)() * spatial_metric(ti::i, ti::j)));
}

template <typename DataType, size_t Dim, typename Frame>
tnsr::ii<DataType, Dim, Frame> a_tilde(
    const Scalar<DataType>& conformal_factor_squared,
    const tnsr::ii<DataType, Dim, Frame>& spatial_metric,
    const tnsr::ii<DataType, Dim, Frame>& extrinsic_curvature,
    const Scalar<DataType>& trace_extrinsic_curvature) {
  tnsr::ii<DataType, Dim, Frame> result{};
  Scalar<DataType> buffer{};
  a_tilde(make_not_null(&result), make_not_null(&buffer),
          conformal_factor_squared, spatial_metric, extrinsic_curvature,
          trace_extrinsic_curvature);
  return result;
}
}  // namespace Ccz4

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define FRAME(data) BOOST_PP_TUPLE_ELEM(1, data)
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(2, data)

#define INSTANTIATE(_, data)                                               \
  template void Ccz4::a_tilde(                                             \
      const gsl::not_null<tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>*>  \
          result,                                                          \
      const gsl::not_null<Scalar<DTYPE(data)>*> buffer,                    \
      const Scalar<DTYPE(data)>& conformal_factor_squared,                 \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>& spatial_metric, \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>&                 \
          extrinsic_curvature,                                             \
      const Scalar<DTYPE(data)>& trace_extrinsic_curvature);               \
  template tnsr::ii<DTYPE(data), DIM(data), FRAME(data)> Ccz4::a_tilde(    \
      const Scalar<DTYPE(data)>& conformal_factor_squared,                 \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>& spatial_metric, \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>&                 \
          extrinsic_curvature,                                             \
      const Scalar<DTYPE(data)>& trace_extrinsic_curvature);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3), (Frame::Grid, Frame::Inertial),
                        (double, DataVector))

#undef INSTANTIATE
#undef DTYPE
#undef FRAME
#undef DIM
