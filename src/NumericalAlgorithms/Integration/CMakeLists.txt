# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY Integration)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  GslQuadAdaptive.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  GslQuadAdaptive.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  ErrorHandling
  GSL::gsl
  Utilities
  )

