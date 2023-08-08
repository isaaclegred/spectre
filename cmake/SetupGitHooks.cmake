# Distributed under the MIT License.
# See LICENSE.txt for details.

# Allow disabling the Git hooks because if you are running the code in a
# container the host may not have all the right things installed.
option(
  USE_GIT_HOOKS
  "Set up the git hooks for sanity checks."
  ON)

if(USE_GIT_HOOKS)
  # Check that the source dir is writable. If it is we set up git hooks, if not
  # then there probably won't be any commits anyway...
  EXECUTE_PROCESS(COMMAND test -w ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE CHECK_SOURCE_DIR_WRITABLE_RESULT)

  # The logic is inverted because shell
  if(NOT CHECK_SOURCE_DIR_WRITABLE_RESULT AND EXISTS ${CMAKE_SOURCE_DIR}/.git)
    find_package(ClangFormat)
    find_package(Git)

    # We use several client-side git hooks to ensure commits are correct as
    # early as possible.
    configure_file(
      ${CMAKE_SOURCE_DIR}/tools/Hooks/pre-commit.sh
      ${CMAKE_SOURCE_DIR}/.git/hooks/pre-commit
      @ONLY
      )

    configure_file(
      ${CMAKE_SOURCE_DIR}/tools/Hooks/CheckFileSize.py
      ${CMAKE_SOURCE_DIR}/.git/hooks/CheckFileSize.py
      @ONLY
      )
  endif()
endif()
