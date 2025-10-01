include(FetchContent)

function(find_or_fetch DEPENDENCY)
    if(DEFINED ENV{CODE})
        set(${DEPENDENCY}_SOURCE_DIR $ENV{CODE}/${DEPENDENCY})
        message(STATUS "Using ${DEPENDENCY} in $ENV{CODE}/${DEPENDENCY}")
    else()
        set(REPO_URL https://github.com/itpplasma/${DEPENDENCY}.git)
        set(REMOTE_BRANCH "")

        if(${DEPENDENCY} STREQUAL "libneo")
            if(DEFINED LIBNEO_BRANCH AND NOT "${LIBNEO_BRANCH}" STREQUAL "")
                string(STRIP "${LIBNEO_BRANCH}" REMOTE_BRANCH)
            elseif(DEFINED ENV{LIBNEO_BRANCH} AND NOT "$ENV{LIBNEO_BRANCH}" STREQUAL "")
                string(STRIP "$ENV{LIBNEO_BRANCH}" REMOTE_BRANCH)
            endif()
        endif()

        if(NOT REMOTE_BRANCH STREQUAL "")
            execute_process(
                COMMAND git ls-remote --heads ${REPO_URL} ${REMOTE_BRANCH}
                OUTPUT_VARIABLE _branch_exists
                RESULT_VARIABLE _ls_result
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            if(NOT _ls_result EQUAL 0 OR _branch_exists STREQUAL "")
                message(WARNING "Requested branch ${REMOTE_BRANCH} not found in ${REPO_URL}; falling back to auto-detected branch")
                set(REMOTE_BRANCH "")
            endif()
        endif()

        if(REMOTE_BRANCH STREQUAL "")
            get_branch_or_main(${REPO_URL} REMOTE_BRANCH)
        endif()
        message(STATUS "Using ${DEPENDENCY} branch ${REMOTE_BRANCH} from ${REPO_URL}")

        FetchContent_Declare(
            ${DEPENDENCY}
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            GIT_REPOSITORY ${REPO_URL}
            GIT_TAG ${REMOTE_BRANCH}
        )
        FetchContent_Populate(${DEPENDENCY})
    endif()

    add_subdirectory(${${DEPENDENCY}_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/${DEPENDENCY}
        EXCLUDE_FROM_ALL
    )
endfunction()


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
