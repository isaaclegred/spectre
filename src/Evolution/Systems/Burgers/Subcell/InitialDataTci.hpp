// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DgSubcell/RdmpTciData.hpp"
#include "Evolution/DgSubcell/Tags/ActiveGrid.hpp"
#include "Evolution/DgSubcell/Tags/DataForRdmpTci.hpp"
#include "Evolution/DgSubcell/Tags/Mesh.hpp"
#include "Evolution/Systems/Burgers/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
template <size_t Dim>
class Mesh;
template <typename TagsList>
class Variables;
/// \endcond

namespace Burgers::subcell {
/// \brief Sets the initial RDMP data.
///
/// Used on the subcells after the TCI marked the DG solution as inadmissible.
struct SetInitialRdmpData {
  using argument_tags =
      tmpl::list<Burgers::Tags::U, evolution::dg::subcell::Tags::ActiveGrid,
                 ::domain::Tags::Mesh<1>,
                 evolution::dg::subcell::Tags::Mesh<1>>;
  using return_tags = tmpl::list<evolution::dg::subcell::Tags::DataForRdmpTci>;

  static void apply(
      gsl::not_null<evolution::dg::subcell::RdmpTciData*> rdmp_tci_data,
      const Scalar<DataVector>& u,
      evolution::dg::subcell::ActiveGrid active_grid, const Mesh<1>& dg_mesh,
      const Mesh<1>& subcell_mesh);
};
}  // namespace Burgers::subcell
