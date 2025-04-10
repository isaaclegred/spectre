# Distributed under the MIT License.
# See LICENSE.txt for details.
set(LIBRARY BnsInitialDataAnalyticData)

add_spectre_library(${LIBRARY})

 if (TARGET SpEC::Exporter)
    target_link_libraries(
      ${LIBRARY}
      PUBLIC
      DataStructures
      ErrorHandling
      Exporter
      GeneralRelativity
      Options
      Utilities
    )
   spectre_target_sources(
     ${LIBRARY}
     PRIVATE
     SpecData.cpp
     SpectreData.cpp
     )
   spectre_target_headers(
     ${LIBRARY}
     INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
     HEADERS
     SpecData.hpp
     SpectreData.hpp
     )
   target_link_libraries(
     ${LIBRARY}
     PUBLIC
     GeneralRelativityAnalyticData
     SpEC::Exporter
     )
else()
    spectre_target_sources(
      ${LIBRARY}
      PRIVATE
      SpectreData.cpp
      )
    spectre_target_headers(
      ${LIBRARY}
      INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
      HEADERS
      SpectreData.hpp
      )
    target_link_libraries(
      ${LIBRARY}
      PUBLIC
      DataStructures
      ErrorHandling
      Exporter
      GeneralRelativity
      GeneralRelativityAnalyticData
      Options
      Utilities
    )
 endif()
