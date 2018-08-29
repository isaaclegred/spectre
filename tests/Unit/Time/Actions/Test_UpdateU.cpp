// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "tests/Unit/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "Time/Actions/UpdateU.hpp"
// IWYU pragma: no_include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags.hpp"
#include "Time/Time.hpp"
#include "Time/TimeSteppers/RungeKutta3.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "tests/Unit/ActionTesting.hpp"

namespace {
struct Var : db::SimpleTag {
  static std::string name() noexcept { return "Var"; }
  using type = double;
};

struct System {
  using variables_tag = Var;
};

struct Metavariables;
struct component
    : ActionTesting::MockArrayComponent<Metavariables, int,
                                        tmpl::list<CacheTags::TimeStepper>> {
  using initial_databox = db::DataBox<tmpl::list<>>;
};

struct Metavariables {
  using system = System;
  using component_list = tmpl::list<component>;
  using const_global_cache_tag_list = tmpl::list<>;
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Time.Actions.UpdateU", "[Unit][Time][Actions]") {
  ActionTesting::ActionRunner<Metavariables> runner{
    {std::make_unique<TimeSteppers::RungeKutta3>()}};
  using variables_tag = Var;
  using dt_variables_tag = Tags::dt<Var>;

  const Slab slab(1., 3.);
  const TimeDelta time_step = slab.duration() / 2;

  using history_tag =
      Tags::HistoryEvolvedVariables<variables_tag, dt_variables_tag>;

  const auto rhs =
      [](const double t, const double y) { return 2. * t - 2. * (y - t * t); };

  auto box =
      db::create<db::AddSimpleTags<Tags::TimeStep, variables_tag, history_tag>>(
          time_step, 1., history_tag::type{});

  const std::array<Time, 3> substep_times{
    {slab.start(), slab.start() + time_step, slab.start() + time_step / 2}};
  // The exact answer is y = x^2, but the integrator would need a
  // smaller step size to get that accurately.
  const std::array<double, 3> expected_values{{3., 3., 10./3.}};

  for (size_t substep = 0; substep < 3; ++substep) {
    db::mutate<history_tag>(
        make_not_null(&box),
        [&rhs, &substep, &substep_times ](
            const gsl::not_null<db::item_type<history_tag>*> history,
            const double& vars) noexcept {
          const Time& time = gsl::at(substep_times, substep);
          history->insert(time, vars, rhs(time.value(), vars));
        },
        db::get<variables_tag>(box));

    box = std::get<0>(runner.apply<component, Actions::UpdateU>(box, 0));

    CHECK(db::get<variables_tag>(box) ==
          approx(gsl::at(expected_values, substep)));
  }
}
