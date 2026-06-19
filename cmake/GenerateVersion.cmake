# Regenerate version.f90 at BUILD time from `git describe`, not just at configure
# time. Run as a script (cmake -P) from a custom command that always executes, so
# the baked-in version follows the current HEAD and dirty state instead of going
# stale after a commit without a reconfigure.
#
# Inputs (passed with -D): SRC (version.f90.in), DST (version.f90), GIT_DIR
# (the source tree to describe).

execute_process(
    COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${GIT_DIR}
    OUTPUT_VARIABLE SIMPLE_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT SIMPLE_VERSION)
    set(SIMPLE_VERSION "unknown")
endif()

# Only rewrite when the version actually changed, so an unchanged HEAD does not
# touch version.f90 and force a needless recompile every build.
set(NEW_CONTENT "")
if(EXISTS ${DST})
    file(READ ${DST} EXISTING_CONTENT)
else()
    set(EXISTING_CONTENT "")
endif()
configure_file(${SRC} ${DST}.tmp @ONLY)
file(READ ${DST}.tmp NEW_CONTENT)
if(NOT "${NEW_CONTENT}" STREQUAL "${EXISTING_CONTENT}")
    file(RENAME ${DST}.tmp ${DST})
    message(STATUS "SIMPLE version: ${SIMPLE_VERSION}")
else()
    file(REMOVE ${DST}.tmp)
endif()
