include(FetchContent)

# Resolve a first-party itpplasma dependency and add it to the build.
#
# Ref precedence (first match wins):
#   1. <DEP>_REF      cache or env: a branch, tag or commit, validated against the
#                     remote and ignored if absent. An upstream release sets this
#                     to build this code against a candidate ref.
#   2. <DEP>_GIT_TAG  cache: an explicit pinned tag (e.g. reference builds).
#   3. <DEP>_RELEASE  cache: the release branch this code tracks by default.
#   4. the current branch if it exists in the remote, otherwise main.
#
# The dependency is always fetched at the resolved ref; there is no local-path
# shortcut.
function(find_or_fetch DEPENDENCY)
    set(REPO_URL https://github.com/itpplasma/${DEPENDENCY}.git)
    string(TOUPPER ${DEPENDENCY} _DEP)

    set(_ref "")
    set(_override "")
    if(DEFINED ${_DEP}_REF AND NOT "${${_DEP}_REF}" STREQUAL "")
        set(_override "${${_DEP}_REF}")
    elseif(DEFINED ENV{${_DEP}_REF} AND NOT "$ENV{${_DEP}_REF}" STREQUAL "")
        set(_override "$ENV{${_DEP}_REF}")
    endif()
    if(NOT "${_override}" STREQUAL "")
        execute_process(
            COMMAND git ls-remote ${REPO_URL} ${_override}
            OUTPUT_VARIABLE _found
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(NOT "${_found}" STREQUAL "")
            set(_ref "${_override}")
        else()
            message(WARNING "${_DEP}_REF='${_override}' not found in ${REPO_URL}; ignoring")
        endif()
    endif()
    if("${_ref}" STREQUAL "" AND DEFINED ${_DEP}_GIT_TAG)
        set(_ref "${${_DEP}_GIT_TAG}")
    endif()
    if("${_ref}" STREQUAL "" AND DEFINED ${_DEP}_RELEASE AND NOT "${${_DEP}_RELEASE}" STREQUAL "")
        set(_ref "${${_DEP}_RELEASE}")
    endif()
    if("${_ref}" STREQUAL "")
        get_branch_or_main(${REPO_URL} _ref)
    endif()
    message(STATUS "Using ${DEPENDENCY} ref ${_ref} from ${REPO_URL}")

    # Fetched first-party dependencies are linked as libraries, not test hosts.
    set(LIBNEO_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(LIBNEO_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(LIBNEO_ENABLE_GOLDEN_TESTS OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        ${DEPENDENCY}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        GIT_REPOSITORY ${REPO_URL}
        GIT_TAG ${_ref}
    )
    FetchContent_GetProperties(${DEPENDENCY})
    if(NOT ${DEPENDENCY}_POPULATED)
        FetchContent_Populate(${DEPENDENCY})
        add_subdirectory(${${DEPENDENCY}_SOURCE_DIR}
            ${CMAKE_CURRENT_BINARY_DIR}/${DEPENDENCY}
            EXCLUDE_FROM_ALL
        )
    endif()
endfunction()

function(get_branch_or_main REPO_URL OUT)
    if(DEFINED ENV{GITHUB_HEAD_REF} AND NOT "$ENV{GITHUB_HEAD_REF}" STREQUAL "")
        set(_branch "$ENV{GITHUB_HEAD_REF}")
    elseif(DEFINED ENV{GITHUB_REF_NAME} AND NOT "$ENV{GITHUB_REF_NAME}" STREQUAL "")
        set(_branch "$ENV{GITHUB_REF_NAME}")
    else()
        execute_process(
            COMMAND git rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE _branch
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()
    if("${_branch}" STREQUAL "" OR "${_branch}" STREQUAL "HEAD")
        set(${OUT} "main" PARENT_SCOPE)
        return()
    endif()
    execute_process(
        COMMAND git ls-remote --heads ${REPO_URL} ${_branch}
        OUTPUT_VARIABLE _exists
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT "${_exists}" STREQUAL "")
        set(${OUT} "${_branch}" PARENT_SCOPE)
    else()
        set(${OUT} "main" PARENT_SCOPE)
    endif()
endfunction()
