# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBS_TO_LINK
  ConstitutiveRelations
  Convergence
  CoordinateMaps
  DiscontinuousGalerkin
  DomainCreators
  Elasticity
  ElasticityAnalyticData
  ElasticityPointwiseFunctions
  ElasticitySolutions
  Elliptic
  EllipticDg
  Events
  Informer
  IO
  LinearOperators
  Options
  Parallel
  ParallelLinearSolver
  Utilities
  )

function(add_elasticity_executable DIM)
  add_spectre_parallel_executable(
    "SolveElasticity${DIM}D"
    SolveElasticity
    Elliptic/Executables/Elasticity
    "Metavariables<${DIM}>"
    "${LIBS_TO_LINK}"
    )
endfunction(add_elasticity_executable)

add_elasticity_executable(2)
add_elasticity_executable(3)
