// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>

#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "Time/Tags/IsUsingTimeSteppingErrorControl.hpp"

SPECTRE_TEST_CASE("Unit.Time.Tags.IsUsingTimeSteppingErrorControl",
                  "[Unit][Time]") {
  TestHelpers::db::test_simple_tag<Tags::IsUsingTimeSteppingErrorControl>(
      "IsUsingTimeSteppingErrorControl");
}
