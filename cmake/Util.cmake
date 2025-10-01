include(FetchContent)

function(get_branch_or_main REPO_URL REMOTE_BRANCH)
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    set(CANDIDATE_BRANCH "${BRANCH}")

    if(CANDIDATE_BRANCH STREQUAL "HEAD" OR CANDIDATE_BRANCH STREQUAL "")
        if(DEFINED ENV{GITHUB_HEAD_REF} AND NOT "$ENV{GITHUB_HEAD_REF}" STREQUAL "")
            set(CANDIDATE_BRANCH "$ENV{GITHUB_HEAD_REF}")
        elseif(DEFINED ENV{GITHUB_REF_NAME} AND NOT "$ENV{GITHUB_REF_NAME}" STREQUAL "")
            set(CANDIDATE_BRANCH "$ENV{GITHUB_REF_NAME}")
        endif()
    endif()

    if(CANDIDATE_BRANCH STREQUAL "HEAD" OR CANDIDATE_BRANCH STREQUAL "")
        set(CANDIDATE_BRANCH "main")
    endif()

    string(STRIP "${CANDIDATE_BRANCH}" CANDIDATE_BRANCH)

    execute_process(
        COMMAND git ls-remote --heads ${REPO_URL} ${CANDIDATE_BRANCH}
        OUTPUT_VARIABLE BRANCH_EXISTS
        ERROR_VARIABLE GIT_ERROR
        RESULT_VARIABLE LS_REMOTE_RESULT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(LS_REMOTE_RESULT EQUAL 0 AND BRANCH_EXISTS)
        set(${REMOTE_BRANCH} ${CANDIDATE_BRANCH} PARENT_SCOPE)
    else()
        if(NOT CANDIDATE_BRANCH STREQUAL "main")
            message(WARNING "Requested branch ${CANDIDATE_BRANCH} not found in ${REPO_URL}; falling back to main")
        endif()
        set(${REMOTE_BRANCH} "main" PARENT_SCOPE)
    endif()
endfunction()

function(find_or_fetch DEPENDENCY)
    string(TOUPPER "${DEPENDENCY}" _dep_upper)
    string(TOLOWER "${DEPENDENCY}" _dep_lower)

    set(_source_dir "")
    set(_use_local FALSE)

    # Allow explicit override via <DEPENDENCY>_SOURCE_DIR cache variable
    set(_override_var "${_dep_upper}_SOURCE_DIR")
    if(DEFINED ${_override_var})
        set(_candidate "${${_override_var}}")
        if(IS_DIRECTORY "${_candidate}")
            set(_source_dir "${_candidate}")
            set(_use_local TRUE)
        else()
            message(FATAL_ERROR "${_override_var}='${_candidate}' is not a directory")
        endif()
    endif()

    if(NOT _use_local AND DEFINED ENV{CODE})
        set(_candidate "$ENV{CODE}/${_dep_lower}")
        if(IS_DIRECTORY "${_candidate}")
            set(_source_dir "${_candidate}")
            set(_use_local TRUE)
        endif()
    endif()

    if(_use_local)
        if("${_dep_lower}" STREQUAL "libneo")
            set(LIBNEO_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
        endif()
        add_subdirectory("${_source_dir}" "${CMAKE_CURRENT_BINARY_DIR}/${_dep_lower}" EXCLUDE_FROM_ALL)
        return()
    endif()

    if("${_dep_lower}" STREQUAL "libneo")
        set(_repo_url "https://github.com/itpplasma/libneo.git")
        set(_branch "")
        if(DEFINED LIBNEO_BRANCH AND NOT "${LIBNEO_BRANCH}" STREQUAL "")
            set(_branch "${LIBNEO_BRANCH}")
        elseif(DEFINED ENV{LIBNEO_BRANCH} AND NOT "$ENV{LIBNEO_BRANCH}" STREQUAL "")
            set(_branch "$ENV{LIBNEO_BRANCH}")
        endif()
        if(_branch STREQUAL "")
            get_branch_or_main(${_repo_url} _branch)
        endif()

        set(LIBNEO_ENABLE_TESTS OFF CACHE BOOL "" FORCE)

        FetchContent_Declare(${_dep_lower}
            GIT_REPOSITORY ${_repo_url}
            GIT_TAG ${_branch}
            GIT_PROGRESS TRUE
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        )
        FetchContent_MakeAvailable(${_dep_lower})
        return()
    endif()

    message(FATAL_ERROR "find_or_fetch does not know how to fetch dependency '${DEPENDENCY}'")
endfunction()
