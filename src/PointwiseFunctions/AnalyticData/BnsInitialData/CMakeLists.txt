# Distributed under the MIT License.
# See LICENSE.txt for details.
set(LIBRARY BnsInitialDataAnalyticData)

if (TARGET SpEC::Exporter)
  add_spectre_library(${LIBRARY})
else()
  add_spectre_library(${LIBRARY}
    INTERFACE)
endif()


 if (TARGET SpEC::Exporter)
    target_link_libraries(
      ${LIBRARY}
      PUBLIC
      DataStructures
      ErrorHandling
      GeneralRelativity
      Options
      Utilities
    )
   spectre_target_sources(
     ${LIBRARY}
     PRIVATE
     SpecData.cpp
     )
   spectre_target_headers(
     ${LIBRARY}
     INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
     HEADERS
     SpecData.hpp
     )
   target_link_libraries(
     ${LIBRARY}
     PUBLIC
     GeneralRelativityAnalyticData
     SpEC::Exporter
     )
else()
    target_link_libraries(
      ${LIBRARY}
      INTERFACE
      DataStructures
      ErrorHandling
      GeneralRelativity
      GeneralRelativityAnalyticData
      Options
      Utilities
    )
 endif()
