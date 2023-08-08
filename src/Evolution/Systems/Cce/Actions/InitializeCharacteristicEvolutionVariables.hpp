// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/Cce/OptionTags.hpp"
#include "NumericalAlgorithms/Spectral/SwshInterpolation.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Cce {
/// \brief The set of actions for use in the CCE evolution system
namespace Actions {

namespace detail {
CREATE_HAS_TYPE_ALIAS(compute_tags)
CREATE_HAS_TYPE_ALIAS_V(compute_tags)
CREATE_GET_TYPE_ALIAS_OR_DEFAULT(compute_tags)
}  // namespace detail

/*!
 * \ingroup ActionsGroup
 * \brief Initializes the main data storage for the  `CharacteristicEvolution`
 * component, which is the singleton that handles the main evolution system for
 * CCE computations.
 *
 * \details Sets up the \ref DataBoxGroup to be ready to take data from the
 * worldtube component, calculate initial data, and start the hypersurface
 * computations.
 *
 * \ref DataBoxGroup changes:
 * - Modifies: nothing
 * - Adds:
 *  - `metavariables::evolved_coordinates_variables_tag`
 *  -
 * ```
 * db::add_tag_prefix<Tags::dt,
 * metavariables::evolved_coordinates_variables_tag>
 * ```
 *  - `Tags::Variables<metavariables::cce_angular_coordinate_tags>`
 *  - `Tags::Variables<metavariables::cce_scri_tags>`
 *  -
 * ```
 * Tags::Variables<tmpl::append<
 * metavariables::cce_integrand_tags,
 * metavariables::cce_integration_independent_tags,
 * metavariables::cce_temporary_equations_tags>>
 * ```
 *  - `Tags::Variables<metavariables::cce_pre_swsh_derivatives_tags>`
 *  - `Tags::Variables<metavariables::cce_transform_buffer_tags>`
 *  - `Tags::Variables<metavariables::cce_swsh_derivative_tags>`
 *  - `Spectral::Swsh::Tags::SwshInterpolator< Tags::CauchyAngularCoords>`
 *  - `Spectral::Swsh::Tags::SwshInterpolator<Tags::PartiallyFlatAngularCoords>`
 * - Removes: nothing
 */
template <typename Metavariables>
struct InitializeCharacteristicEvolutionVariables {
  using const_global_cache_tags =
      tmpl::list<Tags::LMax, Tags::NumberOfRadialPoints>;

  using boundary_value_variables_tag = ::Tags::Variables<
      tmpl::append<typename Metavariables::cce_boundary_communication_tags,
                   typename Metavariables::cce_gauge_boundary_tags>>;
  using scri_variables_tag =
      ::Tags::Variables<typename Metavariables::cce_scri_tags>;
  using volume_variables_tag = ::Tags::Variables<
      tmpl::append<typename Metavariables::cce_integrand_tags,
                   typename Metavariables::cce_integration_independent_tags,
                   typename Metavariables::cce_temporary_equations_tags>>;
  using pre_swsh_derivatives_variables_tag =
      ::Tags::Variables<typename Metavariables::cce_pre_swsh_derivatives_tags>;
  using transform_buffer_variables_tag =
      ::Tags::Variables<typename Metavariables::cce_transform_buffer_tags>;
  using swsh_derivative_variables_tag =
      ::Tags::Variables<typename Metavariables::cce_swsh_derivative_tags>;
  using angular_coordinates_variables_tag =
      ::Tags::Variables<typename Metavariables::cce_angular_coordinate_tags>;
  using coordinate_variables_tag =
      typename Metavariables::evolved_coordinates_variables_tag;
  using dt_coordinate_variables_tag =
      db::add_tag_prefix<::Tags::dt, coordinate_variables_tag>;
  using evolved_swsh_variables_tag =
      ::Tags::Variables<tmpl::list<typename Metavariables::evolved_swsh_tag>>;
  using evolved_swsh_dt_variables_tag =
      db::add_tag_prefix<::Tags::dt, evolved_swsh_variables_tag>;
  using ccm_tag = ::Tags::Variables<typename Metavariables::ccm_psi0>;

  using simple_tags_for_evolution = tmpl::list<
      boundary_value_variables_tag, coordinate_variables_tag,
      dt_coordinate_variables_tag, evolved_swsh_variables_tag,
      evolved_swsh_dt_variables_tag, angular_coordinates_variables_tag,
      scri_variables_tag, volume_variables_tag,
      pre_swsh_derivatives_variables_tag, transform_buffer_variables_tag,
      swsh_derivative_variables_tag,
      Spectral::Swsh::Tags::SwshInterpolator<Tags::CauchyAngularCoords>,
      Spectral::Swsh::Tags::SwshInterpolator<Tags::PartiallyFlatAngularCoords>,
      ccm_tag>;
  using simple_tags =
      tmpl::append<StepChoosers::step_chooser_simple_tags<Metavariables, true>,
                   simple_tags_for_evolution>;

  using compute_tags = tmpl::remove_duplicates<tmpl::join<
      tmpl::transform<typename Metavariables::cce_step_choosers,
                      tmpl::bind<detail::get_compute_tags_or_default_t,
                                 tmpl::_1, tmpl::pin<tmpl::list<>>>>>>;

  template <typename DbTags, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    initialize_impl(make_not_null(&box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }

  template <typename TagList>
  static void initialize_impl(const gsl::not_null<db::DataBox<TagList>*> box) {
    const size_t l_max = db::get<Spectral::Swsh::Tags::LMaxBase>(*box);
    const size_t number_of_radial_points =
        db::get<Spectral::Swsh::Tags::NumberOfRadialPointsBase>(*box);
    const size_t boundary_size =
        Spectral::Swsh::number_of_swsh_collocation_points(l_max);
    const size_t volume_size = boundary_size * number_of_radial_points;
    const size_t transform_buffer_size =
        number_of_radial_points *
        Spectral::Swsh::size_of_libsharp_coefficient_vector(l_max);
    Initialization::mutate_assign<simple_tags_for_evolution>(
        box, typename boundary_value_variables_tag::type{boundary_size},
        typename coordinate_variables_tag::type{boundary_size},
        typename dt_coordinate_variables_tag::type{boundary_size},
        typename evolved_swsh_variables_tag::type{volume_size},
        typename evolved_swsh_dt_variables_tag::type{volume_size},
        typename angular_coordinates_variables_tag::type{boundary_size},
        typename scri_variables_tag::type{boundary_size},
        typename volume_variables_tag::type{volume_size},
        typename pre_swsh_derivatives_variables_tag::type{volume_size, 0.0},
        typename transform_buffer_variables_tag::type{transform_buffer_size,
                                                      0.0},
        typename swsh_derivative_variables_tag::type{volume_size, 0.0},
        Spectral::Swsh::SwshInterpolator{}, Spectral::Swsh::SwshInterpolator{},
        typename ccm_tag::type{boundary_size});
  }
};

}  // namespace Actions
}  // namespace Cce
