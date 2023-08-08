// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Time/TimeSteppers/AdamsBashforth.hpp"

#include <algorithm>
#include <boost/container/small_vector.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <cstddef>
#include <iterator>
#include <limits>
#include <pup.h>
#include <utility>

#include "NumericalAlgorithms/Interpolation/LagrangePolynomial.hpp"
#include "Time/ApproximateTime.hpp"
#include "Time/BoundaryHistory.hpp"
#include "Time/EvolutionOrdering.hpp"
#include "Time/History.hpp"
#include "Time/SelfStart.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Math.hpp"

namespace TimeSteppers {

namespace {
template <typename Iter>
struct TimeFromRecord {
  Time operator()(typename std::iterator_traits<Iter>::reference record) const {
    return record.time_step_id.step_time();
  }
};

// This must be templated on the iterator type rather than the math
// wrapper type because of quirks in the template deduction rules.
template <typename Iter>
auto history_time_iterator(const Iter& it) {
  return boost::transform_iterator(it, TimeFromRecord<Iter>{});
}

// A vector holding one entry per order of integration.
template <typename T>
using OrderVector =
    boost::container::static_vector<T, AdamsBashforth::maximum_order>;

OrderVector<double> constant_coefficients(const size_t order) {
  switch (order) {
    case 0: return {};
    case 1: return {1.};
    case 2: return {1.5, -0.5};
    case 3: return {23.0 / 12.0, -4.0 / 3.0, 5.0 / 12.0};
    case 4: return {55.0 / 24.0, -59.0 / 24.0, 37.0 / 24.0, -3.0 / 8.0};
    case 5: return {1901.0 / 720.0, -1387.0 / 360.0, 109.0 / 30.0,
          -637.0 / 360.0, 251.0 / 720.0};
    case 6: return {4277.0 / 1440.0, -2641.0 / 480.0, 4991.0 / 720.0,
          -3649.0 / 720.0, 959.0 / 480.0, -95.0 / 288.0};
    case 7: return {198721.0 / 60480.0, -18637.0 / 2520.0, 235183.0 / 20160.0,
          -10754.0 / 945.0, 135713.0 / 20160.0, -5603.0 / 2520.0,
          19087.0 / 60480.0};
    case 8: return {16083.0 / 4480.0, -1152169.0 / 120960.0, 242653.0 / 13440.0,
          -296053.0 / 13440.0, 2102243.0 / 120960.0, -115747.0 / 13440.0,
          32863.0 / 13440.0, -5257.0 / 17280.0};
    default:
      ERROR("Bad order: " << order);
  }
}

// This includes the overall factor of step size.
//
// Only T=double is used, but this can be used with T=Rational to
// generate coefficient tables.
template <typename T>
OrderVector<T> variable_coefficients(const OrderVector<T>& control_times) {
  // The argument vector contains the control times for a step from
  // the last time in the list to t=0.

  // The coefficients are, for each j,
  // \int_0^{step} dt ell_j(t; -control_times),
  // where the step size is step=-control_times.back().

  const size_t order = control_times.size();
  OrderVector<T> result;
  for (size_t j = 0; j < order; ++j) {
    // Calculate coefficients of the Lagrange interpolating polynomials,
    // in the standard a_0 + a_1 t + a_2 t^2 + ... form.
    OrderVector<T> poly(order, 0);

    poly[0] = 1;

    for (size_t m = 0; m < order; ++m) {
      if (m == j) {
        continue;
      }
      const T denom =
          1 / (control_times[order - m - 1] - control_times[order - j - 1]);
      for (size_t i = m < j ? m + 1 : m; i > 0; --i) {
        poly[i] =
            (poly[i - 1] + poly[i] * control_times[order - m - 1]) * denom;
      }
      poly[0] *= control_times[order - m - 1] * denom;
    }

    // Integrate p(t), term by term.  We choose the constant of
    // integration so the indefinite integral is zero at t=0.  We do
    // not adjust the indexing, so the t^n term in the integral is in
    // the (n-1)th entry of the vector (as opposed to the nth entry
    // before integrating).
    for (size_t m = 0; m < order; ++m) {
      poly[m] /= m + 1;
    }
    result.push_back(-control_times.back() *
                     evaluate_polynomial(poly, -control_times.back()));
  }
  return result;
}

// Get coefficients for a time step.  Arguments are an iterator pair
// to past times, oldest to newest, and the time step to take.  The
// returned coefficients include the factor of the step size, so, for
// example, the coefficients for Euler's method would be {size}, not
// {1}.
template <typename Iterator, typename Delta>
OrderVector<double> get_coefficients(const Iterator& times_begin,
                                     const Iterator& times_end,
                                     const Delta& step) {
  bool constant_step_size = true;
  OrderVector<double> control_times;
  for (auto t = times_begin; t != times_end; ++t) {
    // Ideally we would also include the slab size in the scale of the
    // roundoff comparison, but there's no good way to get it here,
    // and it should only matter for slabs near time zero.
    if (constant_step_size and not control_times.empty() and
        not equal_within_roundoff(
            t->value() - control_times.back(), step.value(),
            100.0 * std::numeric_limits<double>::epsilon(), abs(t->value()))) {
      constant_step_size = false;
    }
    control_times.push_back(t->value());
  }
  if (constant_step_size) {
    auto result = constant_coefficients(control_times.size());
    alg::for_each(result, [&](double& coef) { coef *= step.value(); });
    return result;
  }

  const double goal_time = control_times.back() + step.value();
  for (auto& t : control_times) {
    t -= goal_time;
  }

  return variable_coefficients(control_times);
}

template <typename T>
void clean_history(const MutableUntypedHistory<T>& history) {
  ASSERT(history.size() >= history.integration_order(),
         "Insufficient data to take an order-" << history.integration_order()
         << " step.  Have " << history.size() << " times, need "
         << history.integration_order());
  while (history.size() > history.integration_order()) {
    history.pop_front();
  }
  if (history.size() > 1) {
    history.discard_value(history[history.size() - 2].time_step_id);
  }
}
}  // namespace

AdamsBashforth::AdamsBashforth(const size_t order) : order_(order) {
  if (order_ < 1 or order_ > maximum_order) {
    ERROR("The order for Adams-Bashforth Nth order must be 1 <= order <= "
          << maximum_order);
  }
}

size_t AdamsBashforth::order() const { return order_; }

size_t AdamsBashforth::error_estimate_order() const { return order_ - 1; }

size_t AdamsBashforth::number_of_past_steps() const { return order_ - 1; }

double AdamsBashforth::stable_step() const {
  if (order_ == 1) {
    return 1.;
  }

  // This is the condition that the characteristic polynomial of the
  // recurrence relation defined by the method has the correct sign at
  // -1.  It is not clear whether this is actually sufficient.
  const auto& coefficients = constant_coefficients(order_);
  double invstep = 0.;
  double sign = 1.;
  for (const auto coef : coefficients) {
    invstep += sign * coef;
    sign = -sign;
  }
  return 1. / invstep;
}

TimeStepId AdamsBashforth::next_time_id(const TimeStepId& current_id,
                                        const TimeDelta& time_step) const {
  ASSERT(current_id.substep() == 0, "Adams-Bashforth should not have substeps");
  return current_id.next_step(time_step);
}

void AdamsBashforth::pup(PUP::er& p) {
  LtsTimeStepper::pup(p);
  p | order_;
}

template <typename T>
void AdamsBashforth::update_u_impl(
    const gsl::not_null<T*> u, const MutableUntypedHistory<T>& history,
    const TimeDelta& time_step) const {
  clean_history(history);
  update_u_common(u, history, time_step, history.integration_order());
}

template <typename T>
bool AdamsBashforth::update_u_impl(
    const gsl::not_null<T*> u, const gsl::not_null<T*> u_error,
    const MutableUntypedHistory<T>& history, const TimeDelta& time_step) const {
  clean_history(history);
  update_u_common(u, history, time_step, history.integration_order());
  // the error estimate is only useful once the history has enough elements to
  // do more than one order of step
  update_u_common(u_error, history, time_step, history.integration_order() - 1);
  *u_error = *u - *u_error;
  return true;
}

template <typename T>
bool AdamsBashforth::dense_update_u_impl(const gsl::not_null<T*> u,
                                         const ConstUntypedHistory<T>& history,
                                         const double time) const {
  const ApproximateTimeDelta time_step{
      time - history.back().time_step_id.step_time().value()};
  update_u_common(u, history, time_step, history.integration_order());
  return true;
}

template <typename T, typename Delta>
void AdamsBashforth::update_u_common(const gsl::not_null<T*> u,
                                     const ConstUntypedHistory<T>& history,
                                     const Delta& time_step,
                                     const size_t order) const {
  ASSERT(
      history.size() > 0,
      "Cannot meaningfully update the evolved variables with an empty history");
  ASSERT(order <= order_,
         "Requested integration order higher than integrator order");

  const auto history_start =
      history.end() -
      static_cast<typename ConstUntypedHistory<T>::difference_type>(order);
  const auto coefficients =
      get_coefficients(history_time_iterator(history_start),
                       history_time_iterator(history.end()), time_step);

  *u = *history.back().value;
  auto coefficient = coefficients.rbegin();
  for (auto history_entry = history_start;
       history_entry != history.end();
       ++history_entry, ++coefficient) {
    *u += *coefficient * history_entry->derivative;
  }
}

template <typename T>
bool AdamsBashforth::can_change_step_size_impl(
    const TimeStepId& time_id, const ConstUntypedHistory<T>& history) const {
  // We need to forbid local time-stepping before initialization is
  // complete.  The self-start procedure itself should never consider
  // changing the step size, but we need to wait during the main
  // evolution until the self-start history has been replaced with
  // "real" values.
  const evolution_less<Time> less{time_id.time_runs_forward()};
  return not ::SelfStart::is_self_starting(time_id) and
         (history.size() == 0 or
          (less(history.back().time_step_id.step_time(),
                time_id.step_time()) and
           std::is_sorted(history_time_iterator(history.begin()),
                          history_time_iterator(history.end()), less)));
}

template <typename T>
void AdamsBashforth::add_boundary_delta_impl(
    const gsl::not_null<T*> result,
    const TimeSteppers::BoundaryHistoryEvaluator<T>& coupling,
    const TimeSteppers::BoundaryHistoryCleaner& cleaner,
    const TimeDelta& time_step) const {
  const auto signed_order =
      static_cast<typename decltype(cleaner.local_end())::difference_type>(
          cleaner.integration_order());

  ASSERT(cleaner.local_size() >= cleaner.integration_order(),
         "Insufficient data to take an order-" << cleaner.integration_order()
         << " step.  Have " << cleaner.local_size() << " times, need "
         << cleaner.integration_order());
  cleaner.local_mark_unneeded(cleaner.local_end() - signed_order);

  if (std::equal(cleaner.local_begin(), cleaner.local_end(),
                 cleaner.remote_end() - signed_order)) {
    // GTS
    ASSERT(cleaner.remote_size() >= cleaner.integration_order(),
           "Insufficient data to take an order-" << cleaner.integration_order()
           << " step.  Have " << cleaner.remote_size() << " times, need "
           << cleaner.integration_order());
    cleaner.remote_mark_unneeded(cleaner.remote_end() - signed_order);
  } else {
    const auto remote_step_for_step_start =
        std::upper_bound(cleaner.remote_begin(), cleaner.remote_end(),
                         *(cleaner.local_end() - 1),
                         evolution_less<Time>{time_step.is_positive()});
    ASSERT(remote_step_for_step_start - cleaner.remote_begin() >= signed_order,
           "Insufficient data to take an order-" << cleaner.integration_order()
           << " step.  Have "
           << remote_step_for_step_start - cleaner.remote_begin()
           << " times before the step, need " << cleaner.integration_order());
    cleaner.remote_mark_unneeded(remote_step_for_step_start - signed_order);
  }

  boundary_impl(result, coupling, *(cleaner.local_end() - 1) + time_step);
}

template <typename T>
void AdamsBashforth::boundary_dense_output_impl(
    const gsl::not_null<T*> result,
    const TimeSteppers::BoundaryHistoryEvaluator<T>& coupling,
    const double time) const {
  if ((coupling.local_end() - 1)->value() == time) {
    // Nothing to do.  The requested time is the start of the step,
    // which is the input value of `result`.
    return;
  }
  return boundary_impl(result, coupling, ApproximateTime{time});
}

namespace {
template <typename T>
class SmallStepIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = Time;
  using pointer = const Time*;
  using reference = const Time&;
  using difference_type = std::ptrdiff_t;

  enum class Side { Local, Remote, Both };

  SmallStepIterator() = default;

  SmallStepIterator(const bool time_runs_forward,
                    typename BoundaryHistoryEvaluator<T>::iterator local_begin,
                    typename BoundaryHistoryEvaluator<T>::iterator remote_begin,
                    typename BoundaryHistoryEvaluator<T>::iterator local_end,
                    typename BoundaryHistoryEvaluator<T>::iterator remote_end)
      : before_{time_runs_forward},
        local_time_(std::move(local_begin)),
        remote_time_(std::move(remote_begin)),
        local_end_(std::move(local_end)),
        remote_end_(std::move(remote_end)) {}

  reference operator*() const {
    return std::max(*local_time_, *remote_time_, before_);
  }
  pointer operator->() const { return &**this; }

  Side side() const {
    if (before_(*local_time_, *remote_time_)) {
      return Side::Remote;
    } else if (before_(*remote_time_, *local_time_)) {
      return Side::Local;
    } else {
      return Side::Both;
    }
  }

  // These are the m^s(n) in the paper.
  const typename BoundaryHistoryEvaluator<T>::iterator& local_iterator() const {
    return local_time_;
  }
  const typename BoundaryHistoryEvaluator<T>::iterator& remote_iterator()
      const {
    return remote_time_;
  }

  SmallStepIterator& operator++() {
    ASSERT(local_time_ != local_end_ and remote_time_ != remote_end_,
           "Overran iterator");
    auto local_candidate = std::next(local_time_);
    auto remote_candidate = std::next(remote_time_);

    // NOLINTNEXTLINE(bugprone-branch-clone)
    if (local_candidate == local_end_ and remote_candidate == remote_end_) {
      local_time_ = std::move(local_candidate);
      remote_time_ = std::move(remote_candidate);
    // NOLINTNEXTLINE(bugprone-branch-clone)
    } else if (local_candidate == local_end_) {
      remote_time_ = std::move(remote_candidate);
    // NOLINTNEXTLINE(bugprone-branch-clone)
    } else if (remote_candidate == remote_end_) {
      local_time_ = std::move(local_candidate);
    } else if (before_(*local_candidate, *remote_candidate)) {
      local_time_ = std::move(local_candidate);
    } else if (before_(*remote_candidate, *local_candidate)) {
      remote_time_ = std::move(remote_candidate);
    } else {
      local_time_ = std::move(local_candidate);
      remote_time_ = std::move(remote_candidate);
    }

    return *this;
  }

  bool done() const {
    return local_time_ == local_end_ and remote_time_ == remote_end_;
  }

 private:
  evolution_less<Time> before_{};
  typename BoundaryHistoryEvaluator<T>::iterator local_time_{};
  typename BoundaryHistoryEvaluator<T>::iterator remote_time_{};
  typename BoundaryHistoryEvaluator<T>::iterator local_end_{};
  typename BoundaryHistoryEvaluator<T>::iterator remote_end_{};
};

template <typename T>
bool operator==(const SmallStepIterator<T>& a, const SmallStepIterator<T>& b) {
  if (a.done() and b.done()) {
    return true;
  }
  if (a.done() or b.done()) {
    return false;
  }
  return *a == *b;
}

template <typename T>
bool operator!=(const SmallStepIterator<T>& a, const SmallStepIterator<T>& b) {
  return not(a == b);
}

template <typename T>
bool operator<(const SmallStepIterator<T>& a, const SmallStepIterator<T>& b) {
  return a.local_iterator() < b.local_iterator() or
         a.remote_iterator() < b.remote_iterator();
}

template <typename T>
bool operator>(const SmallStepIterator<T>& a, const SmallStepIterator<T>& b) {
  return b < a;
}

template <typename It>
It bounded_next_impl(std::input_iterator_tag /*meta*/, It it, const It& bound,
                     typename std::iterator_traits<It>::difference_type n) {
  ASSERT(n >= 0, "Can't advance an input iterator backwards.");
  for (typename std::iterator_traits<It>::difference_type i = 0;
       i < n and it != bound; ++i, ++it) {
  }
  return it;
}

template <typename It>
It bounded_next_impl(std::random_access_iterator_tag /*meta*/, const It& it,
                     typename std::iterator_traits<It>::difference_type n,
                     const It& bound) {
  if (bound - it < n) {
    return bound;
  } else {
    return it + n;
  }
}

template <typename It>
It bounded_next(const It& it, const It& bound, const size_t n) {
  return bounded_next_impl(
      typename std::iterator_traits<It>::iterator_category{}, it, bound,
      static_cast<typename std::iterator_traits<It>::difference_type>(n));
}
}  // namespace

template <typename T, typename TimeType>
void AdamsBashforth::boundary_impl(const gsl::not_null<T*> result,
                                   const BoundaryHistoryEvaluator<T>& coupling,
                                   const TimeType& end_time) const {
  // Might be different from order_ during self-start.
  const auto current_order = coupling.integration_order();

  ASSERT(current_order <= order_,
         "Local history is too long for target order (" << current_order
         << " should not exceed " << order_ << ")");
  ASSERT(coupling.remote_size() >= current_order,
         "Remote history is too short (" << coupling.remote_size()
         << " should be at least " << current_order << ")");

  // Avoid billions of casts
  const auto order_s = static_cast<
      typename BoundaryHistoryEvaluator<T>::iterator::difference_type>(
      current_order);

  // Start and end of the step we are trying to take
  const Time start_time = *(coupling.local_end() - 1);
  const auto time_step = end_time - start_time;

  // We define the local_begin and remote_begin variables as the start
  // of the part of the history relevant to this calculation.
  // Boundary history cleanup happens immediately before the step, but
  // boundary dense output happens before that, so there may be data
  // left over that was needed for the previous step and has not been
  // cleaned out yet.
  const auto local_begin = coupling.local_end() - order_s;

  if (std::equal(local_begin, coupling.local_end(),
                 coupling.remote_end() - order_s)) {
    // No local time-stepping going on.
    const auto coefficients =
        get_coefficients(local_begin, coupling.local_end(), time_step);

    auto local_it = local_begin;
    auto remote_it = coupling.remote_end() - order_s;
    for (auto coefficients_it = coefficients.rbegin();
         coefficients_it != coefficients.rend();
         ++coefficients_it, ++local_it, ++remote_it) {
      *result += *coefficients_it * *coupling(local_it, remote_it);
    }
    return;
  }

  ASSERT(current_order == order_,
         "Cannot perform local time-stepping while self-starting.");

  const evolution_less<> less{time_step.is_positive()};
  const auto remote_begin =
      std::upper_bound(coupling.remote_begin(), coupling.remote_end(),
                       start_time, less) -
      order_s;

  ASSERT(std::is_sorted(local_begin, coupling.local_end(), less),
         "Local history not in order");
  ASSERT(std::is_sorted(remote_begin, coupling.remote_end(), less),
         "Remote history not in order");
  ASSERT(not less(start_time, *(remote_begin + (order_s - 1))),
         "Remote history does not extend far enough back");
  ASSERT(less(*(coupling.remote_end() - 1), end_time),
         "Please supply only older data: " << *(coupling.remote_end() - 1)
         << " is not before " << end_time);

  using difference_type = std::ptrdiff_t;

  SmallStepIterator<T> contributing_small_step(
      time_step.is_positive(), local_begin, remote_begin, coupling.local_end(),
      coupling.remote_end());
  SmallStepIterator<T> small_step_of_current_step = std::next(
      contributing_small_step, static_cast<difference_type>(current_order - 1));
  while (*small_step_of_current_step != start_time) {
    ++contributing_small_step;
    ++small_step_of_current_step;
  }

  // The size of this vector is the number of small steps in the
  // current step, which will almost always be 1 or 2, but could be
  // larger if two elements decide to do 4:1 stepping or something.
  boost::container::small_vector<OrderVector<double>, 2>
      small_step_coefficients{};
  {
    auto coefficient_eval_begin = contributing_small_step;
    auto coefficient_eval_end = small_step_of_current_step;
    while (coefficient_eval_end != SmallStepIterator<T>{}) {
      auto next_end = std::next(coefficient_eval_end);
      const double next_time = next_end == SmallStepIterator<T>{}
                                   ? end_time.value()
                                   : next_end->value();
      small_step_coefficients.push_back(get_coefficients(
          coefficient_eval_begin, next_end,
          ApproximateTimeDelta{next_time - coefficient_eval_end->value()}));
      ++coefficient_eval_begin;
      coefficient_eval_end = std::move(next_end);
    }
  }

  // Sum over the small steps that contribute to this step, doing the
  // appropriate interpolation for each.
  size_t current_step_minus_contributing_step = current_order - 1;
  for (; contributing_small_step != SmallStepIterator<T>{};
       ++contributing_small_step, --current_step_minus_contributing_step) {
    if (contributing_small_step.side() != SmallStepIterator<T>::Side::Local) {
      double overall_prefactor = 0.0;
      auto small_step_within_current_step = static_cast<size_t>(std::max(
          -static_cast<difference_type>(current_step_minus_contributing_step),
          static_cast<difference_type>(0)));
      const size_t small_step_within_current_step_end =
          std::min(small_step_coefficients.size(),
                   current_order - current_step_minus_contributing_step);
      for (;
           small_step_within_current_step < small_step_within_current_step_end;
           ++small_step_within_current_step) {
        overall_prefactor +=
            small_step_coefficients[small_step_within_current_step]
                                   [small_step_within_current_step +
                                    current_step_minus_contributing_step];
      }
      if (contributing_small_step.side() == SmallStepIterator<T>::Side::Both) {
        *result += overall_prefactor *
                   *coupling(contributing_small_step.local_iterator(),
                             contributing_small_step.remote_iterator());
      } else {
        // Side::Remote
        OrderVector<double> past_steps(current_order);
        std::transform(
            coupling.local_end() - static_cast<difference_type>(current_order),
            coupling.local_end(), past_steps.begin(),
            [](const Time& t) { return t.value(); });
        size_t interpolation_index = 0;
        for (auto interpolation_time =
                 coupling.local_end() -
                 static_cast<difference_type>(current_order);
             interpolation_time != coupling.local_end();
             ++interpolation_time, ++interpolation_index) {
          const double coefficient =
              overall_prefactor *
              lagrange_polynomial(interpolation_index,
                                  contributing_small_step->value(),
                                  past_steps.begin(), past_steps.end());
          *result += coefficient *
                     *coupling(interpolation_time,
                               contributing_small_step.remote_iterator());
        }
      }
    } else {
      // Side::Local
      auto interpolation_time =
          std::max(small_step_of_current_step.remote_iterator(),
                   contributing_small_step.remote_iterator()) -
          static_cast<difference_type>(current_order - 1);
      const auto interpolation_time_end =
          bounded_next(contributing_small_step, SmallStepIterator<T>{},
                       current_order)
              .remote_iterator();
      for (; interpolation_time != interpolation_time_end;
           ++interpolation_time) {
        double coefficient = 0.0;
        size_t small_step_within_current_step_index = 0;
        auto small_step_within_current_step = small_step_of_current_step;
        while (small_step_within_current_step.remote_iterator() <
               interpolation_time) {
          ++small_step_within_current_step;
          ++small_step_within_current_step_index;
        }
        auto small_step_within_current_step_end = contributing_small_step;
        const auto bound_from_interpolation_time = bounded_next(
            interpolation_time, coupling.remote_end(), current_order);
        for (size_t i = 0;
             i < current_order and
             small_step_within_current_step_end != SmallStepIterator<T>{} and
             small_step_within_current_step_end.remote_iterator() <
                 bound_from_interpolation_time;
             ++i) {
          ++small_step_within_current_step_end;
        }
        auto remote_steps_end = coupling.remote_begin();
        double lagrange_factor = std::numeric_limits<double>::signaling_NaN();
        for (; small_step_within_current_step !=
               small_step_within_current_step_end;
             ++small_step_within_current_step,
             ++small_step_within_current_step_index) {
          if (remote_steps_end !=
              small_step_within_current_step.remote_iterator() + 1) {
            remote_steps_end =
                small_step_within_current_step.remote_iterator() + 1;
            OrderVector<double> past_steps(current_order);
            std::transform(
                remote_steps_end - static_cast<difference_type>(current_order),
                remote_steps_end, past_steps.begin(),
                [](const Time& t) { return t.value(); });
            lagrange_factor = lagrange_polynomial(
                current_order -
                    static_cast<size_t>(remote_steps_end - interpolation_time),
                contributing_small_step->value(), past_steps.begin(),
                past_steps.end());
          }
          coefficient +=
              lagrange_factor *
              small_step_coefficients[small_step_within_current_step_index]
                                     [small_step_within_current_step_index +
                                      current_step_minus_contributing_step];
        }
        *result +=
            coefficient * *coupling(contributing_small_step.local_iterator(),
                                    interpolation_time);
      }
    }
  }
}

bool operator==(const AdamsBashforth& lhs, const AdamsBashforth& rhs) {
  return lhs.order_ == rhs.order_;
}

bool operator!=(const AdamsBashforth& lhs, const AdamsBashforth& rhs) {
  return not(lhs == rhs);
}

TIME_STEPPER_DEFINE_OVERLOADS(AdamsBashforth)
LTS_TIME_STEPPER_DEFINE_OVERLOADS(AdamsBashforth)
}  // namespace TimeSteppers

PUP::able::PUP_ID TimeSteppers::AdamsBashforth::my_PUP_ID = 0;  // NOLINT
