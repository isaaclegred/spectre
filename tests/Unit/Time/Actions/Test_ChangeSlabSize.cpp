// Distributed under the MIT License.
// See LICENSE.txt for details.

// The event portion of slab size changing relies on a reduction, so
// it cannot currently be tested with the mocking code.  This tests
// the action portion.

#include "Framework/TestingFramework.hpp"

#include <cstdint>
#include <initializer_list>
#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Framework/ActionTesting.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"  // IWYU pragma: keep
#include "Time/Actions/ChangeSlabSize.hpp"
#include "Time/AdaptiveSteppingDiagnostics.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags.hpp"
#include "Time/Tags/AdaptiveSteppingDiagnostics.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/Rk3HesthavenSsp.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// IWYU pragma: no_forward_declare ActionTesting::InitializeDataBox
namespace Tags {
template <typename Tag>
struct Next;
template <typename Tag>
struct dt;
}  // namespace Tags

namespace {
struct Var : db::SimpleTag {
  using type = double;
};

struct Component;

struct Metavariables {
  using component_list = tmpl::list<Component>;
};

struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;

  using const_global_cache_tags = tmpl::list<Tags::TimeStepper<TimeStepper>>;

  using simple_tags =
      tmpl::list<Tags::TimeStepId, Tags::Next<Tags::TimeStepId>, Tags::TimeStep,
                 Tags::Next<Tags::TimeStep>, Tags::HistoryEvolvedVariables<Var>,
                 Tags::AdaptiveSteppingDiagnostics>;
  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<
                     Parallel::Phase::Initialization,
                     tmpl::list<ActionTesting::InitializeDataBox<simple_tags>>>,
                 Parallel::PhaseActions<Parallel::Phase::Testing,
                                        tmpl::list<Actions::ChangeSlabSize>>>;
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Time.Actions.ChangeSlabSize", "[Unit][Time][Actions]") {
  register_classes_with_charm<TimeSteppers::Rk3HesthavenSsp>();

  ActionTesting::MockRuntimeSystem<Metavariables> runner{
      {std::make_unique<TimeSteppers::Rk3HesthavenSsp>()}};

  ActionTesting::emplace_component_and_initialize<Component>(&runner, 0, {});
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

  auto& box = ActionTesting::get_databox<Component>(make_not_null(&runner), 0);

  for (const bool time_runs_forward : {true, false}) {
    Slab slab(1.5, 2.0);
    Time start_time;

    const auto resize_slab = [&slab, &start_time,
                              &time_runs_forward](const double length) {
      if (time_runs_forward) {
        slab = slab.with_duration_from_start(length);
        start_time = slab.start();
      } else {
        slab = slab.with_duration_to_end(length);
        start_time = slab.end();
      }
    };

    resize_slab(slab.duration().value());

    const auto get_step = [](const TimeStepId& id) {
      return (id.time_runs_forward() ? 1 : -1) *
             id.step_time().slab().duration() / 2;
    };

    db::mutate<Tags::TimeStepId, Tags::Next<Tags::TimeStepId>, Tags::TimeStep,
               Tags::Next<Tags::TimeStep>, Tags::AdaptiveSteppingDiagnostics>(
        make_not_null(&box),
        [&get_step, &start_time, &time_runs_forward](
            const gsl::not_null<TimeStepId*> id,
            const gsl::not_null<TimeStepId*> next_id,
            const gsl::not_null<TimeDelta*> step,
            const gsl::not_null<TimeDelta*> next_step,
            const gsl::not_null<AdaptiveSteppingDiagnostics*> diags,
            const TimeStepper& stepper) {
          *id = TimeStepId(time_runs_forward, 3, start_time);
          *step = get_step(*id);
          *next_step = get_step(*id);
          *next_id = stepper.next_time_id(*id, *step);
          *diags = AdaptiveSteppingDiagnostics{1, 2, 3, 4, 5};
        },
        db::get<Tags::TimeStepper<>>(box));

    using ExpectedMessages =
        ChangeSlabSize_detail::NumberOfExpectedMessagesInbox;
    using NewSize = ChangeSlabSize_detail::NewSlabSizeInbox;
    auto& inboxes = runner.inboxes<Component>().at(0);

    const auto check_box = [&box, &get_step](const TimeStepId& id,
                                             const uint64_t changes) {
      CHECK(db::get<Tags::TimeStepId>(box) == id);
      CHECK(db::get<Tags::TimeStep>(box) == get_step(id));
      CHECK(db::get<Tags::Next<Tags::TimeStepId>>(box) ==
            db::get<Tags::TimeStepper<>>(box).next_time_id(
                db::get<Tags::TimeStepId>(box), db::get<Tags::TimeStep>(box)));
      CHECK(db::get<Tags::AdaptiveSteppingDiagnostics>(box) ==
            AdaptiveSteppingDiagnostics{1, 2 + changes, 3, 4, 5});
    };

    // Nothing to do
    {
      runner.next_action<Component>(0);
      check_box(TimeStepId(time_runs_forward, 3, start_time), 0);
    }

    // Simple case
    {
      get<ExpectedMessages>(inboxes)[3].insert(ExpectedMessages::NoData{});
      REQUIRE_FALSE(ActionTesting::next_action_if_ready<Component>(
          make_not_null(&runner), 0));
      get<NewSize>(inboxes)[3].insert(1.0);
      runner.next_action<Component>(0);
      resize_slab(1.0);
      check_box(TimeStepId(time_runs_forward, 3, start_time), 1);
      CHECK(get<ExpectedMessages>(inboxes).empty());
      CHECK(get<NewSize>(inboxes).empty());
      get<ExpectedMessages>(inboxes).clear();
      get<NewSize>(inboxes).clear();
    }

    // Multiple messages at multiple times
    {
      get<ExpectedMessages>(inboxes)[3].insert(ExpectedMessages::NoData{});
      get<ExpectedMessages>(inboxes)[3].insert(ExpectedMessages::NoData{});
      get<ExpectedMessages>(inboxes)[4].insert(ExpectedMessages::NoData{});
      REQUIRE_FALSE(ActionTesting::next_action_if_ready<Component>(
          make_not_null(&runner), 0));
      get<NewSize>(inboxes)[3].insert(2.0);
      REQUIRE_FALSE(ActionTesting::next_action_if_ready<Component>(
          make_not_null(&runner), 0));
      get<NewSize>(inboxes)[4].insert(0.5);
      REQUIRE_FALSE(ActionTesting::next_action_if_ready<Component>(
          make_not_null(&runner), 0));
      get<ExpectedMessages>(inboxes)[4].insert(ExpectedMessages::NoData{});
      get<NewSize>(inboxes)[3].insert(3.0);
      runner.next_action<Component>(0);
      resize_slab(2.0);
      check_box(TimeStepId(time_runs_forward, 3, start_time), 2);
      CHECK(get<ExpectedMessages>(inboxes).size() == 1);
      CHECK(get<ExpectedMessages>(inboxes).count(4) == 1);
      CHECK(get<NewSize>(inboxes).size() == 1);
      CHECK(get<NewSize>(inboxes).count(4) == 1);
      get<ExpectedMessages>(inboxes).clear();
      get<NewSize>(inboxes).clear();
    }

    // Check interior of slab
    {
      db::mutate<Tags::TimeStepId, Tags::Next<Tags::TimeStepId>>(
          make_not_null(&box),
          [&start_time, &time_runs_forward](
              const gsl::not_null<TimeStepId*> id,
              const gsl::not_null<TimeStepId*> next_id, const TimeDelta& step,
              const TimeStepper& stepper) {
            *id = TimeStepId(time_runs_forward, 3, start_time + step);
            *next_id = stepper.next_time_id(*id, step);
          },
          db::get<Tags::TimeStep>(box), db::get<Tags::TimeStepper<>>(box));
      const TimeStepId initial_id = db::get<Tags::TimeStepId>(box);
      get<ExpectedMessages>(inboxes)[3].insert(ExpectedMessages::NoData{});
      get<ExpectedMessages>(inboxes)[4].insert(ExpectedMessages::NoData{});
      runner.next_action<Component>(0);
      check_box(initial_id, 2);
      CHECK(get<ExpectedMessages>(inboxes).size() == 2);
      CHECK(get<ExpectedMessages>(inboxes).count(3) == 1);
      CHECK(get<ExpectedMessages>(inboxes).count(4) == 1);
      CHECK(get<NewSize>(inboxes).empty());
      get<NewSize>(inboxes)[3].insert(0.1);
      get<NewSize>(inboxes)[4].insert(0.1);
      runner.next_action<Component>(0);
      check_box(initial_id, 2);
      get<ExpectedMessages>(inboxes).clear();
      get<NewSize>(inboxes).clear();
    }

    // Check at a substep
    {
      db::mutate<Tags::TimeStepId, Tags::Next<Tags::TimeStepId>>(
          make_not_null(&box),
          [](const gsl::not_null<TimeStepId*> id,
             const gsl::not_null<TimeStepId*> next_id, const TimeDelta& step,
             const TimeStepper& stepper) {
            const auto local_slab = id->step_time().slab();
            while (id->substep_time() != local_slab.start().value() and
                   id->substep_time() != local_slab.end().value()) {
              *id = *next_id;
              REQUIRE(id->substep() != 0);
              *next_id = stepper.next_time_id(*id, step);
            }
          },
          db::get<Tags::TimeStep>(box), db::get<Tags::TimeStepper<>>(box));
      const TimeStepId initial_id = db::get<Tags::TimeStepId>(box);
      get<ExpectedMessages>(inboxes)[3].insert(ExpectedMessages::NoData{});
      get<ExpectedMessages>(inboxes)[4].insert(ExpectedMessages::NoData{});
      runner.next_action<Component>(0);
      check_box(initial_id, 2);
      CHECK(get<ExpectedMessages>(inboxes).size() == 2);
      CHECK(get<ExpectedMessages>(inboxes).count(3) == 1);
      CHECK(get<ExpectedMessages>(inboxes).count(4) == 1);
      CHECK(get<NewSize>(inboxes).empty());
      get<NewSize>(inboxes)[3].insert(0.1);
      get<NewSize>(inboxes)[4].insert(0.1);
      runner.next_action<Component>(0);
      check_box(initial_id, 2);
      get<ExpectedMessages>(inboxes).clear();
      get<NewSize>(inboxes).clear();
    }
  }
}
