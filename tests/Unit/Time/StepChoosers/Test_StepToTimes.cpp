// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <pup.h>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/Tags/Metavariables.hpp"
#include "Time/Slab.hpp"
#include "Time/StepChoosers/StepChooser.hpp"
#include "Time/StepChoosers/StepToTimes.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeSequence.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct Metavariables {
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<StepChooser<StepChooserUse::Slab>,
                             tmpl::list<StepChoosers::StepToTimes>>,
                  tmpl::pair<TimeSequence<double>,
                             TimeSequences::all_time_sequences<double>>>;
  };
  using component_list = tmpl::list<>;
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Time.StepChoosers.StepToTimes", "[Unit][Time]") {
  register_factory_classes_with_charm<Metavariables>();

  const auto requested = [](const double now, std::vector<double> times,
                            const double step) {
    double result = -1.0;
    const auto impl = [&result, &step](const TimeStepId& now_id,
                                       const std::vector<double>& impl_times) {
      CAPTURE(now_id);
      CAPTURE(impl_times);
      CAPTURE(step);

      using Specified = TimeSequences::Specified<double>;
      const StepChoosers::StepToTimes step_to_times(
          std::make_unique<Specified>(impl_times));
      const std::unique_ptr<StepChooser<StepChooserUse::Slab>>
          step_to_times_base = std::make_unique<StepChoosers::StepToTimes>(
              std::make_unique<Specified>(impl_times));

      auto box = db::create<db::AddSimpleTags<
          Parallel::Tags::MetavariablesImpl<Metavariables>, Tags::TimeStepId>>(
          Metavariables{}, now_id);

      const auto answer = step_to_times(now_id, step);
      if (result == -1.0) {
        result = answer.first;
      } else {
        CHECK(result == answer.first);
      }
      CHECK(step_to_times_base->desired_step(step, box) ==
            std::make_pair(result, true));
      CHECK(serialize_and_deserialize(step_to_times)(now_id, step) ==
            std::make_pair(result, true));
      CHECK(serialize_and_deserialize(step_to_times_base)
                ->desired_step(step, box) == std::make_pair(result, true));
    };
    impl(TimeStepId(true, 0, Slab(now, now + 1.0).start()), times);
    alg::for_each(times, [](double& x) { return x = -x; });
    impl(TimeStepId(false, 0, Slab(-now, -now + 1.0).start()), times);
    return result;
  };

  static constexpr double infinity = std::numeric_limits<double>::infinity();

  CHECK(requested(3.0, {}, infinity) == infinity);
  CHECK(requested(3.0, {1.0}, infinity) == infinity);
  CHECK(requested(3.0, {3.0}, infinity) == infinity);
  CHECK(requested(3.0, {5.0}, infinity) == 2.0);
  {
    const double request = requested(3.0, {5.0}, 1.0);
    CHECK(request > 1.0);
    CHECK(request < 2.0);
  }
  {
    const double request = requested(3.0, {5.0}, 1.9375);
    CHECK(request > 1.0);
    CHECK(request < 1.9375);
  }
  {
    const double request = requested(3.0, {5.0}, 1.0625);
    CHECK(request > 1.0);
    CHECK(request < 2.0);
  }

  const auto one_five_check = [&requested](const std::vector<double>& times) {
    CHECK(requested(0.0, times, infinity) == 1.0);
    CHECK(requested(1.0, times, infinity) == 4.0);
    CHECK(requested(2.0, times, infinity) == 3.0);
    CHECK(requested(5.0, times, infinity) == infinity);
    CHECK(requested(7.0, times, infinity) == infinity);
  };

  one_five_check({1.0, 5.0});
  one_five_check({5.0, 1.0});
  one_five_check({1.0, 5.0, 1.0, 5.0, 1.0, 5.0});

  const auto check_rounding = [&requested](const double start,
                                           const double step) {
    CAPTURE(start);
    CAPTURE(step);
    auto scaled_approx = Approx::custom().scale(std::abs(start));
    CHECK(requested(std::nextafter(start, +infinity), {start, start + step},
                    infinity) == scaled_approx(step));
    CHECK(requested(std::nextafter(start, -infinity), {start, start + step},
                    infinity) == scaled_approx(step));
  };

  check_rounding(0.0, 1.0);
  check_rounding(1.0, 1.0);
  check_rounding(-1.0, 1.0);
  check_rounding(1.0e5, 1.0);
  check_rounding(-1.0e5, 1.0);
  check_rounding(1.0, 1.0e5);
  check_rounding(-1.0, 1.0e5);
  check_rounding(1.0e5, 1.0e5);
  check_rounding(-1.0e5, 1.0e5);

  TestHelpers::test_creation<std::unique_ptr<StepChooser<StepChooserUse::Slab>>,
                             Metavariables>(
      "StepToTimes:\n"
      "  Times:\n"
      "    Specified:\n"
      "      Values: [5.0, 3.0, 6.0]");

  CHECK(not StepChoosers::StepToTimes{}.uses_local_data());
}
