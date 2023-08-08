// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTriggers/Times.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeSequence.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct Metavariables {
  using component_list = tmpl::list<>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<DenseTrigger, tmpl::list<DenseTriggers::Times>>,
                  tmpl::pair<TimeSequence<double>,
                             TimeSequences::all_time_sequences<double>>>;
  };
};

void check_one_direction(const std::vector<double>& trigger_times,
                         const double current_time,
                         const bool expected_is_triggered,
                         const double expected_next_check,
                         const bool time_runs_forward) {
  CAPTURE(time_runs_forward);

  std::stringstream creation_string;
  creation_string.precision(std::numeric_limits<double>::max_digits10);
  creation_string << "Times:\n"
                  << "  Specified:\n"
                  << "    Values: [";
  for (auto time : trigger_times) {
    creation_string << time << ",";
  }
  creation_string << "]";
  CAPTURE(creation_string.str());

  const auto trigger = serialize_and_deserialize(
      TestHelpers::test_creation<std::unique_ptr<DenseTrigger>, Metavariables>(
          creation_string.str()));

  const Slab slab(1.0e10, 2.0e10);
  auto box = db::create<
      db::AddSimpleTags<Parallel::Tags::MetavariablesImpl<Metavariables>,
                        Tags::TimeStepId, Tags::Time>>(
      Metavariables{}, TimeStepId(time_runs_forward, 100, slab.start()),
      current_time);
  Parallel::GlobalCache<Metavariables> cache{};
  const int array_index = 0;
  const void* component = nullptr;

  CHECK_FALSE(trigger->previous_trigger_time().has_value());
  CHECK(trigger->is_triggered(box, cache, array_index, component) ==
        expected_is_triggered);
  CHECK(trigger->next_check_time(box, cache, array_index, component) ==
        std::optional{expected_next_check});
  if (expected_is_triggered == std::optional{true}) {
    CHECK_FALSE(trigger->previous_trigger_time().has_value());
    db::mutate<::Tags::Time>(
        [](const gsl::not_null<double*> time) { *time += 0.01; },
        make_not_null(&box));
    CHECK(trigger->is_triggered(box, cache, array_index, component) ==
          std::optional{false});
    REQUIRE(trigger->previous_trigger_time().has_value());
    CHECK(trigger->previous_trigger_time().value() == current_time);
  } else {
    CHECK_FALSE(trigger->previous_trigger_time().has_value());
    db::mutate<::Tags::Time>(
        [](const gsl::not_null<double*> time) { *time += 0.01; },
        make_not_null(&box));
    CHECK(trigger->is_triggered(box, cache, array_index, component) ==
          std::optional{false});
    CHECK_FALSE(trigger->previous_trigger_time().has_value());
  }
}

void check_both_directions(std::vector<double> trigger_times,
                           const double current_time,
                           const bool expected_is_triggered,
                           const double expected_next_check) {
  check_one_direction(trigger_times, current_time, expected_is_triggered,
                      expected_next_check, true);
  alg::for_each(trigger_times, [](double& t) { t = -t; });
  check_one_direction(trigger_times, -current_time, expected_is_triggered,
                      -expected_next_check, false);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.EventsAndDenseTriggers.DenseTriggers.Times",
                  "[Unit][Evolution]") {
  register_factory_classes_with_charm<Metavariables>();

  const auto infinity = std::numeric_limits<double>::infinity();

  check_both_directions({}, 1.0, false, infinity);

  check_both_directions({1.7}, 1.6, false, 1.7);
  check_both_directions({1.7}, 1.7, true, infinity);
  check_both_directions({1.7}, 1.8, false, infinity);

  check_both_directions({1.7, 1.9}, 1.6, false, 1.7);
  check_both_directions({1.7, 1.9}, 1.7, true, 1.9);
  check_both_directions({1.7, 1.9}, 1.8, false, 1.9);
  check_both_directions({1.7, 1.9}, 1.9, true, infinity);
  check_both_directions({1.7, 1.9}, 2.0, false, infinity);
}
