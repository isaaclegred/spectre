// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Utilities/NoSuchType.hpp"
#include "Utilities/TMPL.hpp"

// Check if SpEC is linked and therefore we can load SpEC  data

#ifdef HAS_SPEC_EXPORTER
#include "PointwiseFunctions/AnalyticData/BnsInitialData/SpecData.hpp"
using SpecDataList = tmpl::list<BnsInitialData::AnalyticData::SpecData<1>>;
#else
using SpecDataList = NoSuchType;
#endif

namespace BnsInitialData::AnalyticData {

using all_analytic_data =
    tmpl::conditional_t<std::is_same_v<SpecDataList, NoSuchType>, tmpl::list<>,
                        SpecDataList>;
}  // namespace BnsInitialData::AnalyticData
