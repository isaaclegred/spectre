// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/ProjectionOperators.hpp"

#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace gr {
template <typename DataType, size_t VolumeDim, typename Frame>
tnsr::II<DataType, VolumeDim, Frame> transverse_projection_operator(
    const tnsr::II<DataType, VolumeDim, Frame>& inverse_spatial_metric,
    const tnsr::I<DataType, VolumeDim, Frame>& normal_vector) {
  tnsr::II<DataType, VolumeDim, Frame> projection_tensor(
      get_size(get<0>(normal_vector)));
  transverse_projection_operator(make_not_null(&projection_tensor),
                                 inverse_spatial_metric, normal_vector);
  return projection_tensor;
}

template <typename DataType, size_t VolumeDim, typename Frame>
void transverse_projection_operator(
    const gsl::not_null<tnsr::II<DataType, VolumeDim, Frame>*>
        projection_tensor,
    const tnsr::II<DataType, VolumeDim, Frame>& inverse_spatial_metric,
    const tnsr::I<DataType, VolumeDim, Frame>& normal_vector) {
  for (size_t i = 0; i < VolumeDim; ++i) {
    for (size_t j = i; j < VolumeDim; ++j) {
      projection_tensor->get(i, j) =
          inverse_spatial_metric.get(i, j) -
          normal_vector.get(i) * normal_vector.get(j);
    }
  }
}

template <typename DataType, size_t VolumeDim, typename Frame>
tnsr::ii<DataType, VolumeDim, Frame> transverse_projection_operator(
    const tnsr::ii<DataType, VolumeDim, Frame>& spatial_metric,
    const tnsr::i<DataType, VolumeDim, Frame>& normal_one_form) {
  tnsr::ii<DataType, VolumeDim, Frame> projection_tensor(
      get_size(get<0>(normal_one_form)));
  transverse_projection_operator(make_not_null(&projection_tensor),
                                 spatial_metric, normal_one_form);
  return projection_tensor;
}

template <typename DataType, size_t VolumeDim, typename Frame>
void transverse_projection_operator(
    const gsl::not_null<tnsr::ii<DataType, VolumeDim, Frame>*>
        projection_tensor,
    const tnsr::ii<DataType, VolumeDim, Frame>& spatial_metric,
    const tnsr::i<DataType, VolumeDim, Frame>& normal_one_form) {
  for (size_t i = 0; i < VolumeDim; ++i) {
    for (size_t j = i; j < VolumeDim; ++j) {
      projection_tensor->get(i, j) =
          spatial_metric.get(i, j) -
          normal_one_form.get(i) * normal_one_form.get(j);
    }
  }
}

template <typename DataType, size_t VolumeDim, typename Frame>
tnsr::Ij<DataType, VolumeDim, Frame> transverse_projection_operator(
    const tnsr::I<DataType, VolumeDim, Frame>& normal_vector,
    const tnsr::i<DataType, VolumeDim, Frame>& normal_one_form) {
  tnsr::Ij<DataType, VolumeDim, Frame> projection_tensor(
      get_size(get<0>(normal_vector)));
  transverse_projection_operator(make_not_null(&projection_tensor),
                                 normal_vector, normal_one_form);
  return projection_tensor;
}

template <typename DataType, size_t VolumeDim, typename Frame>
void transverse_projection_operator(
    const gsl::not_null<tnsr::Ij<DataType, VolumeDim, Frame>*>
        projection_tensor,
    const tnsr::I<DataType, VolumeDim, Frame>& normal_vector,
    const tnsr::i<DataType, VolumeDim, Frame>& normal_one_form) {
  for (size_t i = 0; i < VolumeDim; ++i) {
    for (size_t j = 0; j < VolumeDim; ++j) {
      projection_tensor->get(i, j) =
          -normal_vector.get(i) * normal_one_form.get(j);
    }
    projection_tensor->get(i, i) += 1.;
  }
}

template <typename DataType, size_t VolumeDim, typename Frame>
tnsr::AA<DataType, VolumeDim, Frame> transverse_projection_operator(
    const tnsr::AA<DataType, VolumeDim, Frame>& inverse_spacetime_metric,
    const tnsr::A<DataType, VolumeDim, Frame>& spacetime_normal_vector,
    const tnsr::I<DataType, VolumeDim, Frame>& interface_unit_normal_vector) {
  tnsr::AA<DataType, VolumeDim, Frame> projection_tensor(
      get_size(get<0>(spacetime_normal_vector)));
  transverse_projection_operator(
      make_not_null(&projection_tensor), inverse_spacetime_metric,
      spacetime_normal_vector, interface_unit_normal_vector);
  return projection_tensor;
}

template <typename DataType, size_t VolumeDim, typename Frame>
void transverse_projection_operator(
    const gsl::not_null<tnsr::AA<DataType, VolumeDim, Frame>*>
        projection_tensor,
    const tnsr::AA<DataType, VolumeDim, Frame>& inverse_spacetime_metric,
    const tnsr::A<DataType, VolumeDim, Frame>& spacetime_normal_vector,
    const tnsr::I<DataType, VolumeDim, Frame>& interface_unit_normal_vector) {
  for (size_t a = 0, b = 0; b < VolumeDim + 1; ++b) {
    projection_tensor->get(a, b) =
        inverse_spacetime_metric.get(a, b) +
        spacetime_normal_vector.get(a) * spacetime_normal_vector.get(b);
  }
  for (size_t a = 1; a < VolumeDim + 1; ++a) {
    for (size_t b = a; b < VolumeDim + 1; ++b) {
      projection_tensor->get(a, b) =
          inverse_spacetime_metric.get(a, b) +
          spacetime_normal_vector.get(a) * spacetime_normal_vector.get(b) -
          interface_unit_normal_vector.get(a - 1) *
              interface_unit_normal_vector.get(b - 1);
    }
  }
}

template <typename DataType, size_t VolumeDim, typename Frame>
tnsr::aa<DataType, VolumeDim, Frame> transverse_projection_operator(
    const tnsr::aa<DataType, VolumeDim, Frame>& spacetime_metric,
    const tnsr::a<DataType, VolumeDim, Frame>& spacetime_normal_one_form,
    const tnsr::i<DataType, VolumeDim, Frame>& interface_unit_normal_one_form) {
  tnsr::aa<DataType, VolumeDim, Frame> projection_tensor(
      get_size(get<0>(spacetime_normal_one_form)));
  transverse_projection_operator(make_not_null(&projection_tensor),
                                 spacetime_metric, spacetime_normal_one_form,
                                 interface_unit_normal_one_form);
  return projection_tensor;
}

template <typename DataType, size_t VolumeDim, typename Frame>
void transverse_projection_operator(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame>*>
        projection_tensor,
    const tnsr::aa<DataType, VolumeDim, Frame>& spacetime_metric,
    const tnsr::a<DataType, VolumeDim, Frame>& spacetime_normal_one_form,
    const tnsr::i<DataType, VolumeDim, Frame>& interface_unit_normal_one_form) {
  for (size_t a = 0, b = 0; b < VolumeDim + 1; ++b) {
    projection_tensor->get(a, b) =
        spacetime_metric.get(a, b) +
        spacetime_normal_one_form.get(a) * spacetime_normal_one_form.get(b);
  }
  for (size_t a = 1; a < VolumeDim + 1; ++a) {
    for (size_t b = a; b < VolumeDim + 1; ++b) {
      projection_tensor->get(a, b) =
          spacetime_metric.get(a, b) +
          spacetime_normal_one_form.get(a) * spacetime_normal_one_form.get(b) -
          interface_unit_normal_one_form.get(a - 1) *
              interface_unit_normal_one_form.get(b - 1);
    }
  }
}

template <typename DataType, size_t VolumeDim, typename Frame>
tnsr::Ab<DataType, VolumeDim, Frame> transverse_projection_operator(
    const tnsr::A<DataType, VolumeDim, Frame>& spacetime_normal_vector,
    const tnsr::a<DataType, VolumeDim, Frame>& spacetime_normal_one_form,
    const tnsr::I<DataType, VolumeDim, Frame>& interface_unit_normal_vector,
    const tnsr::i<DataType, VolumeDim, Frame>& interface_unit_normal_one_form) {
  tnsr::Ab<DataType, VolumeDim, Frame> projection_tensor(
      get_size(get<0>(spacetime_normal_vector)));
  transverse_projection_operator(
      make_not_null(&projection_tensor), spacetime_normal_vector,
      spacetime_normal_one_form, interface_unit_normal_vector,
      interface_unit_normal_one_form);
  return projection_tensor;
}

template <typename DataType, size_t VolumeDim, typename Frame>
void transverse_projection_operator(
    const gsl::not_null<tnsr::Ab<DataType, VolumeDim, Frame>*>
        projection_tensor,
    const tnsr::A<DataType, VolumeDim, Frame>& spacetime_normal_vector,
    const tnsr::a<DataType, VolumeDim, Frame>& spacetime_normal_one_form,
    const tnsr::I<DataType, VolumeDim, Frame>& interface_unit_normal_vector,
    const tnsr::i<DataType, VolumeDim, Frame>& interface_unit_normal_one_form) {
  for (size_t a = 0, b = 0; b < VolumeDim + 1; ++b) {
    projection_tensor->get(a, b) =
        spacetime_normal_vector.get(a) * spacetime_normal_one_form.get(b);
    projection_tensor->get(b, a) =
        spacetime_normal_vector.get(b) * spacetime_normal_one_form.get(a);
  }
  get<0, 0>(*projection_tensor) += 1.;

  for (size_t a = 1; a < VolumeDim + 1; ++a) {
    for (size_t b = 1; b < VolumeDim + 1; ++b) {
      projection_tensor->get(a, b) =
          spacetime_normal_vector.get(a) * spacetime_normal_one_form.get(b) -
          interface_unit_normal_vector.get(a - 1) *
              interface_unit_normal_one_form.get(b - 1);
    }
    projection_tensor->get(a, a) += 1.;
  }
}
}  // namespace gr

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)
#define FRAME(data) BOOST_PP_TUPLE_ELEM(2, data)

#define INSTANTIATE(_, data)                                                 \
  template tnsr::II<DTYPE(data), DIM(data), FRAME(data)>                     \
  gr::transverse_projection_operator(                                        \
      const tnsr::II<DTYPE(data), DIM(data), FRAME(data)>&                   \
          inverse_spatial_metric,                                            \
      const tnsr::I<DTYPE(data), DIM(data), FRAME(data)>& normal_vector);    \
  template tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>                     \
  gr::transverse_projection_operator(                                        \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>& spatial_metric,   \
      const tnsr::i<DTYPE(data), DIM(data), FRAME(data)>& normal_one_form);  \
  template tnsr::Ij<DTYPE(data), DIM(data), FRAME(data)>                     \
  gr::transverse_projection_operator(                                        \
      const tnsr::I<DTYPE(data), DIM(data), FRAME(data)>& normal_vector,     \
      const tnsr::i<DTYPE(data), DIM(data), FRAME(data)>& normal_one_form);  \
  template void gr::transverse_projection_operator(                          \
      const gsl::not_null<tnsr::II<DTYPE(data), DIM(data), FRAME(data)>*>    \
          projection_tensor,                                                 \
      const tnsr::II<DTYPE(data), DIM(data), FRAME(data)>&                   \
          inverse_spatial_metric,                                            \
      const tnsr::I<DTYPE(data), DIM(data), FRAME(data)>& normal_vector);    \
  template void gr::transverse_projection_operator(                          \
      const gsl::not_null<tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>*>    \
          projection_tensor,                                                 \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>& spatial_metric,   \
      const tnsr::i<DTYPE(data), DIM(data), FRAME(data)>& normal_one_form);  \
  template void gr::transverse_projection_operator(                          \
      const gsl::not_null<tnsr::Ij<DTYPE(data), DIM(data), FRAME(data)>*>    \
          projection_tensor,                                                 \
      const tnsr::I<DTYPE(data), DIM(data), FRAME(data)>& normal_vector,     \
      const tnsr::i<DTYPE(data), DIM(data), FRAME(data)>& normal_one_form);  \
  template tnsr::AA<DTYPE(data), DIM(data), FRAME(data)>                     \
  gr::transverse_projection_operator(                                        \
      const tnsr::AA<DTYPE(data), DIM(data), FRAME(data)>&                   \
          inverse_spacetime_metric,                                          \
      const tnsr::A<DTYPE(data), DIM(data), FRAME(data)>&                    \
          spacetime_normal_vector,                                           \
      const tnsr::I<DTYPE(data), DIM(data), FRAME(data)>&                    \
          interface_unit_normal_vector);                                     \
  template tnsr::aa<DTYPE(data), DIM(data), FRAME(data)>                     \
  gr::transverse_projection_operator(                                        \
      const tnsr::aa<DTYPE(data), DIM(data), FRAME(data)>& spacetime_metric, \
      const tnsr::a<DTYPE(data), DIM(data), FRAME(data)>&                    \
          spacetime_normal_one_form,                                         \
      const tnsr::i<DTYPE(data), DIM(data), FRAME(data)>&                    \
          interface_unit_normal_one_form);                                   \
  template tnsr::Ab<DTYPE(data), DIM(data), FRAME(data)>                     \
  gr::transverse_projection_operator(                                        \
      const tnsr::A<DTYPE(data), DIM(data), FRAME(data)>&                    \
          spacetime_normal_vector,                                           \
      const tnsr::a<DTYPE(data), DIM(data), FRAME(data)>&                    \
          spacetime_normal_one_form,                                         \
      const tnsr::I<DTYPE(data), DIM(data), FRAME(data)>&                    \
          interface_unit_normal_vector,                                      \
      const tnsr::i<DTYPE(data), DIM(data), FRAME(data)>&                    \
          interface_unit_normal_one_form);                                   \
  template void gr::transverse_projection_operator(                          \
      const gsl::not_null<tnsr::AA<DTYPE(data), DIM(data), FRAME(data)>*>    \
          projection_tensor,                                                 \
      const tnsr::AA<DTYPE(data), DIM(data), FRAME(data)>&                   \
          inverse_spacetime_metric,                                          \
      const tnsr::A<DTYPE(data), DIM(data), FRAME(data)>&                    \
          spacetime_normal_vector,                                           \
      const tnsr::I<DTYPE(data), DIM(data), FRAME(data)>&                    \
          interface_unit_normal_vector);                                     \
  template void gr::transverse_projection_operator(                          \
      const gsl::not_null<tnsr::aa<DTYPE(data), DIM(data), FRAME(data)>*>    \
          projection_tensor,                                                 \
      const tnsr::aa<DTYPE(data), DIM(data), FRAME(data)>& spacetime_metric, \
      const tnsr::a<DTYPE(data), DIM(data), FRAME(data)>&                    \
          spacetime_normal_one_form,                                         \
      const tnsr::i<DTYPE(data), DIM(data), FRAME(data)>&                    \
          interface_unit_normal_one_form);                                   \
  template void gr::transverse_projection_operator(                          \
      const gsl::not_null<tnsr::Ab<DTYPE(data), DIM(data), FRAME(data)>*>    \
          projection_tensor,                                                 \
      const tnsr::A<DTYPE(data), DIM(data), FRAME(data)>&                    \
          spacetime_normal_vector,                                           \
      const tnsr::a<DTYPE(data), DIM(data), FRAME(data)>&                    \
          spacetime_normal_one_form,                                         \
      const tnsr::I<DTYPE(data), DIM(data), FRAME(data)>&                    \
          interface_unit_normal_vector,                                      \
      const tnsr::i<DTYPE(data), DIM(data), FRAME(data)>&                    \
          interface_unit_normal_one_form);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3), (double, DataVector),
                        (Frame::Grid, Frame::Inertial))

#undef DIM
#undef DTYPE
#undef FRAME
#undef INSTANTIATE
