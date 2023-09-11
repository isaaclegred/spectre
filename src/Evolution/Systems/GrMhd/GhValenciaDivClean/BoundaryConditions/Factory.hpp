// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/BoundaryConditions/ConstraintPreservingFreeOutflow.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/BoundaryConditions/DirichletAnalytic.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/BoundaryConditions/DirichletFreeOutflow.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/BoundaryConditions/Factory.hpp"
#include "Utilities/TMPL.hpp"

namespace grmhd::GhValenciaDivClean::BoundaryConditions {
namespace detail {
template <typename ClassList>
struct remove_periodic_conditions {
  using type = tmpl::remove_if<
      ClassList,
      std::is_base_of<domain::BoundaryConditions::MarkAsPeriodic, tmpl::_1>>;
};

template <typename ClassList>
using remove_periodic_conditions_t =
    typename remove_periodic_conditions<ClassList>::type;
}  // namespace detail

// remove the periodic BCs from the creatable classes of the
// individual systems; for the remaining conditions, include a
// `ProductOfConditions` for each pair with compatible `bc_type`s.
/// Typelist of standard BoundaryConditions
using standard_boundary_conditions =
    tmpl::push_back<detail::remove_periodic_conditions_t<
                        typename grmhd::ValenciaDivClean::BoundaryConditions::
                            standard_boundary_conditions>,
                    domain::BoundaryConditions::Periodic<BoundaryCondition>>;

/// Boundary conditions that work with finite difference.
using standard_fd_boundary_conditions = tmpl::list<
    ConstraintPreservingFreeOutflow,
    domain::BoundaryConditions::Periodic<BoundaryCondition>,
    grmhd::GhValenciaDivClean::BoundaryConditions::DirichletAnalytic,
    grmhd::GhValenciaDivClean::BoundaryConditions::DirichletFreeOutflow>;
}  // namespace grmhd::GhValenciaDivClean::BoundaryConditions
