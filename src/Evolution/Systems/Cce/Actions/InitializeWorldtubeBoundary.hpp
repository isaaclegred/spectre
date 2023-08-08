// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/Cce/AnalyticBoundaryDataManager.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/RobinsonTrautman.hpp"
#include "Evolution/Systems/Cce/BoundaryData.hpp"
#include "Evolution/Systems/Cce/InterfaceManagers/GhInterfaceManager.hpp"
#include "Evolution/Systems/Cce/InterfaceManagers/GhLocalTimeStepping.hpp"
#include "Evolution/Systems/Cce/InterfaceManagers/GhLockstep.hpp"
#include "Evolution/Systems/Cce/OptionTags.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"
#include "NumericalAlgorithms/Interpolation/SpanInterpolator.hpp"
#include "NumericalAlgorithms/Spectral/SwshCollocation.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/Invoke.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/TimeSteppers/LtsTimeStepper.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeString.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Cce {
/// \cond
template <typename Metavariables>
struct H5WorldtubeBoundary;
template <typename Metavariables>
struct AnalyticWorldtubeBoundary;
template <typename Metavariables>
struct GhWorldtubeBoundary;
/// \endcond
namespace Actions {

namespace detail {
template <typename Initializer, typename ManagerTags,
          typename BoundaryCommunicationTagsList>
struct InitializeWorldtubeBoundaryBase {
  using simple_tags_from_options = ManagerTags;
  using const_global_cache_tags = tmpl::list<Tags::LMax>;

  using simple_tags =
      tmpl::list<::Tags::Variables<BoundaryCommunicationTagsList>>;

  template <typename DataBoxTagsList, typename... InboxTags,
            typename ArrayIndex, typename Metavariables, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DataBoxTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    if constexpr (std::is_same_v<Tags::AnalyticBoundaryDataManager,
                                 tmpl::front<ManagerTags>>) {
      if (dynamic_cast<const Solutions::RobinsonTrautman*>(
              &(db::get<Tags::AnalyticBoundaryDataManager>(box)
                    .get_generator())) != nullptr) {
        if(db::get<::Tags::TimeStepper<>>(box).number_of_substeps() != 1) {
          ERROR(
              "Do not use RobinsonTrautman analytic solution with a "
              "substep-based timestepper. This is to prevent severe slowdowns "
              "in the current RobinsonTrautman implementation. See the "
              "documentation for the RobinsonTrautman solution for details.");
        }
      }
    }
    const size_t l_max = db::get<Tags::LMax>(box);
    Variables<BoundaryCommunicationTagsList> boundary_variables{
        Spectral::Swsh::number_of_swsh_collocation_points(l_max)};

    Initialization::mutate_assign<simple_tags>(make_not_null(&box),
                                               std::move(boundary_variables));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace detail

/*!
 * \ingroup ActionsGroup
 * \brief Generic action for initializing various worldtube boundary components.
 *
 * \details See specializations of this class for initialization details for
 * individual worldtube components.
 */
template <typename WorldtubeComponent>
struct InitializeWorldtubeBoundary;

/*!
 * \ingroup ActionsGroup
 * \brief Initializes a H5WorldtubeBoundary
 *
 * \details Uses:
 * - initialization tag
 * `Cce::Tags::H5WorldtubeBoundaryDataManager`,
 * - const global cache tag `Cce::Tags::LMax`.
 *
 * Databox changes:
 * - Adds:
 *   - `Cce::Tags::H5WorldtubeBoundaryDataManager`
 *   - `Tags::Variables<typename
 * Metavariables::cce_boundary_communication_tags>`
 * - Removes: nothing
 * - Modifies: nothing
 */
template <typename Metavariables>
struct InitializeWorldtubeBoundary<H5WorldtubeBoundary<Metavariables>>
    : public detail::InitializeWorldtubeBoundaryBase<
          InitializeWorldtubeBoundary<H5WorldtubeBoundary<Metavariables>>,
          tmpl::list<Tags::H5WorldtubeBoundaryDataManager>,
          typename Metavariables::cce_boundary_communication_tags> {
  using base_type = detail::InitializeWorldtubeBoundaryBase<
      InitializeWorldtubeBoundary<H5WorldtubeBoundary<Metavariables>>,
      tmpl::list<Tags::H5WorldtubeBoundaryDataManager>,
      typename Metavariables::cce_boundary_communication_tags>;
  using base_type::apply;
  using typename base_type::simple_tags;
  using const_global_cache_tags =
      tmpl::list<Tags::LMax, Tags::EndTimeFromFile, Tags::StartTimeFromFile>;
  using typename base_type::simple_tags_from_options;
};

/*!
 * \ingroup ActionsGroup
 * \brief Initializes a GhWorldtubeBoundary
 *
 * \details Uses:
 * - initialization tags
 * `Cce::Tags::GhWorldtubeBoundaryDataManager`, `Tags::GhInterfaceManager`
 * - const global cache tags `Tags::LMax`, `Tags::ExtractionRadius`.
 *
 * Databox changes:
 * - Adds:
 *   - `Tags::Variables<typename
 * Metavariables::cce_boundary_communication_tags>`
 * - Removes: nothing
 * - Modifies: nothing
 */
template <typename Metavariables>
struct InitializeWorldtubeBoundary<GhWorldtubeBoundary<Metavariables>>
    : public detail::InitializeWorldtubeBoundaryBase<
          InitializeWorldtubeBoundary<GhWorldtubeBoundary<Metavariables>>,
          tmpl::list<Tags::GhInterfaceManager,
                     Tags::SelfStartGhInterfaceManager>,
          typename Metavariables::cce_boundary_communication_tags> {
  using base_type = detail::InitializeWorldtubeBoundaryBase<
      InitializeWorldtubeBoundary<GhWorldtubeBoundary<Metavariables>>,
      tmpl::list<Tags::GhInterfaceManager, Tags::SelfStartGhInterfaceManager>,
      typename Metavariables::cce_boundary_communication_tags>;
  using base_type::apply;
  using typename base_type::simple_tags;

  using const_global_cache_tags =
      tmpl::list<Tags::LMax, InitializationTags::ExtractionRadius,
                 Tags::NoEndTime, Tags::SpecifiedStartTime>;
  using typename base_type::simple_tags_from_options;
};

/*!
 * \ingroup ActionsGroup
 * \brief Initializes an AnalyticWorldtubeBoundary
 *
 * \details Uses:
 * - initialization tag
 * `Cce::Tags::AnalyticBoundaryDataManager`,
 * - const global cache tags `Tags::LMax`,
 * `Tags::ExtractionRadius`.
 *
 * Databox changes:
 * - Adds:
 *   - `Tags::Variables<typename
 *     Metavariables::cce_boundary_communication_tags>`
 * - Removes: nothing
 * - Modifies: nothing
 */
template <typename Metavariables>
struct InitializeWorldtubeBoundary<AnalyticWorldtubeBoundary<Metavariables>>
    : public detail::InitializeWorldtubeBoundaryBase<
          InitializeWorldtubeBoundary<AnalyticWorldtubeBoundary<Metavariables>>,
          tmpl::list<Tags::AnalyticBoundaryDataManager,
                     Tags::CceEvolutionPrefix<::Tags::TimeStepper<
                         tmpl::conditional_t<Metavariables::local_time_stepping,
                                             LtsTimeStepper, TimeStepper>>>>,
          typename Metavariables::cce_boundary_communication_tags> {
  using base_type = detail::InitializeWorldtubeBoundaryBase<
      InitializeWorldtubeBoundary<AnalyticWorldtubeBoundary<Metavariables>>,
      tmpl::list<Tags::AnalyticBoundaryDataManager,
                 Tags::CceEvolutionPrefix<::Tags::TimeStepper<
                     tmpl::conditional_t<Metavariables::local_time_stepping,
                                         LtsTimeStepper, TimeStepper>>>>,
      typename Metavariables::cce_boundary_communication_tags>;
  using base_type::apply;
  using typename base_type::simple_tags;
  using const_global_cache_tags =
      tmpl::list<Tags::LMax, Tags::SpecifiedEndTime, Tags::SpecifiedStartTime>;
  using typename base_type::simple_tags_from_options;
};
}  // namespace Actions
}  // namespace Cce
