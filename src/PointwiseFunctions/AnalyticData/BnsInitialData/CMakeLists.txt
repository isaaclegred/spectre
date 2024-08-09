# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY BnsInitialDataAnalyticData)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  TovStar.cpp
  SpecInitialData.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  TovStar.hpp
  SpecInitialData.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  ErrorHandling
  GeneralRelativity
  Options
  RelativisticEulerSolutions
  Utilities
  )
