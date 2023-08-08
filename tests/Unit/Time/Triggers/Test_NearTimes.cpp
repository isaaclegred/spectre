// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <memory>
#include <pup.h>
#include <pup_stl.h>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/Tags/Metavariables.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Time.hpp"
#include "Time/TimeSequence.hpp"
#include "Time/Triggers/NearTimes.hpp"
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
        tmpl::map<tmpl::pair<TimeSequence<double>,
                             TimeSequences::all_time_sequences<double>>,
                  tmpl::pair<Trigger, tmpl::list<Triggers::NearTimes>>>;
  };
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Time.Triggers.NearTimes", "[Unit][Time]") {
  register_factory_classes_with_charm<Metavariables>();

  using Direction = Triggers::NearTimes::Direction;
  using Unit = Triggers::NearTimes::Unit;

  const auto check = [](const double time_in, const double range_in,
                        const std::vector<double>& trigger_times_in,
                        const Direction direction, const bool expected) {
    const auto check_signs = [&direction, &expected, &trigger_times_in,
                              &time_in](const TimeDelta time_step_in,
                                        const double range, const Unit unit) {
      const auto check_calls = [&direction, &expected, &range, &unit](
                                   std::vector<double> trigger_times,
                                   const double time,
                                   const TimeDelta& time_step) {
        CAPTURE(trigger_times);
        CAPTURE(range);
        CAPTURE(static_cast<int>(unit));
        CAPTURE(static_cast<int>(direction));
        CAPTURE(time);
        CAPTURE(time_step);
        const std::unique_ptr<Trigger> trigger =
            std::make_unique<Triggers::NearTimes>(
                std::make_unique<TimeSequences::Specified<double>>(
                    std::move(trigger_times)),
                range, unit, direction);
        const auto sent_trigger = serialize_and_deserialize(trigger);

        const auto box = db::create<
            db::AddSimpleTags<Parallel::Tags::MetavariablesImpl<Metavariables>,
                              Tags::Time, Tags::TimeStep>>(
            Metavariables{}, time, time_step);

        CHECK(trigger->is_triggered(box) == expected);
        CHECK(sent_trigger->is_triggered(box) == expected);
      };

      check_calls(trigger_times_in, time_in, time_step_in);
      std::vector<double> negated_trigger_times = trigger_times_in;
      alg::for_each(negated_trigger_times, [](double& x) { x = -x; });
      check_calls(negated_trigger_times, -time_in, -time_step_in);
    };

    // Absolute position of slab should not matter.
    // Slab size 100, step size 10
    const TimeDelta large_step = Slab(10000.0, 10100.0).duration() / 10;
    // Slab size 0.1, step_size 0.1
    const TimeDelta small_step = Slab(10000.0, 10000.1).duration();

    check_signs(large_step, range_in, Unit::Time);
    check_signs(large_step, range_in / 10.0, Unit::Step);
    check_signs(large_step, range_in / 100.0, Unit::Slab);
    check_signs(small_step, range_in, Unit::Time);
    check_signs(small_step, range_in / 0.1, Unit::Step);
    check_signs(small_step, range_in / 0.1, Unit::Slab);
  };

  check(5.0, 3.0, {}, Direction::Both, false);
  check(5.0, 3.0, {}, Direction::Before, false);
  check(5.0, 3.0, {}, Direction::After, false);
  check(5.0, 3.0, {6.0}, Direction::Both, true);
  check(5.0, 3.0, {6.0}, Direction::Before, true);
  check(5.0, 3.0, {6.0}, Direction::After, false);
  check(5.0, 3.0, {4.0}, Direction::Both, true);
  check(5.0, 3.0, {4.0}, Direction::Before, false);
  check(5.0, 3.0, {4.0}, Direction::After, true);
  check(5.0, 3.0, {4.0, 6.0}, Direction::Both, true);
  check(5.0, 3.0, {4.0, 6.0}, Direction::Before, true);
  check(5.0, 3.0, {4.0, 6.0}, Direction::After, true);
  check(5.0, 3.0, {9.0}, Direction::Both, false);
  check(5.0, 3.0, {9.0}, Direction::Before, false);
  check(5.0, 3.0, {9.0}, Direction::After, false);
  check(5.0, 3.0, {1.0}, Direction::Both, false);
  check(5.0, 3.0, {1.0}, Direction::Before, false);
  check(5.0, 3.0, {1.0}, Direction::After, false);
  check(5.0, 3.0, {1.0, 9.0}, Direction::Both, false);
  check(5.0, 3.0, {1.0, 9.0}, Direction::Before, false);
  check(5.0, 3.0, {1.0, 9.0}, Direction::After, false);

  TestHelpers::test_creation<std::unique_ptr<Trigger>, Metavariables>(
      "NearTimes:\n"
      "  Times:\n"
      "    Specified:\n"
      "      Values: [2.0, 1.0, 3.0, 2.0]\n"
      "  Range: 0.3\n"
      "  Direction: Before\n"
      "  Unit: Time");
  CHECK(TestHelpers::test_creation<Unit>("Time") == Unit::Time);
  CHECK(TestHelpers::test_creation<Unit>("Step") == Unit::Step);
  CHECK(TestHelpers::test_creation<Unit>("Slab") == Unit::Slab);
  CHECK(TestHelpers::test_creation<Direction>("Before") == Direction::Before);
  CHECK(TestHelpers::test_creation<Direction>("After") == Direction::After);
  CHECK(TestHelpers::test_creation<Direction>("Both") == Direction::Both);
}
