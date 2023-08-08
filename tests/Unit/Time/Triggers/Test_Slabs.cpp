// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstdint>
#include <initializer_list>
#include <memory>
#include <pup.h>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/Tags/Metavariables.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeSequence.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/Triggers/Slabs.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct Metavariables {
  using component_list = tmpl::list<>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<TimeSequence<std::uint64_t>,
                             TimeSequences::all_time_sequences<std::uint64_t>>,
                  tmpl::pair<Trigger, tmpl::list<Triggers::Slabs>>>;
  };
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Time.Triggers.Slabs", "[Unit][Time]") {
  register_factory_classes_with_charm<Metavariables>();

  const auto trigger =
      TestHelpers::test_creation<std::unique_ptr<Trigger>, Metavariables>(
          "Slabs:\n"
          "  Specified:\n"
          "    Values: [3, 6, 8]");

  const auto sent_trigger = serialize_and_deserialize(trigger);

  const Slab slab(0., 1.);
  auto box = db::create<db::AddSimpleTags<
      Parallel::Tags::MetavariablesImpl<Metavariables>, Tags::TimeStepId>>(
      Metavariables{}, TimeStepId(true, 0, slab.start()));
  for (const bool expected :
       {false, false, false, true, false, false, true, false, true, false}) {
    CHECK(sent_trigger->is_triggered(box) == expected);
    db::mutate<Tags::TimeStepId>(
        make_not_null(&box), [&slab](const gsl::not_null<TimeStepId*> time_id) {
          *time_id = TimeStepId(true, time_id->slab_number(),
                                time_id->step_time(), 1, slab.duration(),
                                time_id->step_time().value());
        });
    CHECK_FALSE(sent_trigger->is_triggered(box));
    db::mutate<Tags::TimeStepId>(
        make_not_null(&box), [](const gsl::not_null<TimeStepId*> time_id) {
          *time_id = TimeStepId(true, time_id->slab_number() + 1,
                                time_id->step_time());
        });
  }
}
