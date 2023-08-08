// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/ForceFree/Fluxes.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ForceFree::subcell {
namespace detail {
/*!
 * \brief Helper function that calls `ForceFree::Fluxes::apply` by
 * retrieving the return and argument tags from the Variables object `vars`.
 */
template <typename TagsList, typename... ReturnTags, typename... ArgumentTags>
void compute_fluxes_impl(const gsl::not_null<Variables<TagsList>*> vars,
                         tmpl::list<ReturnTags...> /*meta*/,
                         tmpl::list<ArgumentTags...> /*meta*/) {
  ForceFree::Fluxes::apply(make_not_null(&get<ReturnTags>(*vars))...,
                           get<ArgumentTags>(*vars)...);
}
}  // namespace detail

/*!
 * \brief Compute fluxes for subcell variables.
 */
template <typename TagsList>
void compute_fluxes(const gsl::not_null<Variables<TagsList>*> vars) {
  detail::compute_fluxes_impl(vars, typename ForceFree::Fluxes::return_tags{},
                              typename ForceFree::Fluxes::argument_tags{});
}
}  // namespace ForceFree::subcell
