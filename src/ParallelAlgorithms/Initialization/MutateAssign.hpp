// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace Initialization {
namespace detail {
template <typename... MutateTags, typename BoxTags, typename... Args>
SPECTRE_ALWAYS_INLINE constexpr void mutate_assign_impl(
    // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
    const gsl::not_null<db::DataBox<BoxTags>*> box,
    tmpl::list<MutateTags...> /*meta*/, Args&&... args) {
  static_assert(sizeof...(MutateTags) == sizeof...(args),
                "The number of arguments passed to `mutate_assign` must be "
                "equal to the number of tags passed.");
  db::mutate<MutateTags...>(
      [&args...](const auto... box_args) {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
        EXPAND_PACK_LEFT_TO_RIGHT((*box_args = std::forward<Args>(args)));
      },
      box);
}
}  // namespace detail

/*!
 * \ingroup DataBoxGroup
 * \brief Perform a mutation to the \ref DataBoxGroup `box`, assigning the
 * `args` to the tags in `MutateTagList` in order.
 */
template <typename MutateTagList, typename BoxTags, typename... Args>
SPECTRE_ALWAYS_INLINE constexpr void mutate_assign(
    // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
    const gsl::not_null<db::DataBox<BoxTags>*> box, Args&&... args) {
  // The impl works for zero tags, but we can skip it to improve
  // compilation performance.
  if constexpr (tmpl::size<MutateTagList>::value != 0) {
    detail::mutate_assign_impl(box, MutateTagList{},
                               std::forward<Args>(args)...);
  } else {
    (void)box;
    expand_pack(args...);
  }
}
}  // namespace Initialization
