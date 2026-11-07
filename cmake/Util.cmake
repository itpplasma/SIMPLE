include(FetchContent)

function(find_or_fetch DEPENDENCY)
    string(TOUPPER "${DEPENDENCY}" _dep_upper)
    string(TOLOWER "${DEPENDENCY}" _dep_lower)

    set(_source_dir "")
    set(_use_local FALSE)

    # Allow explicit override via <DEPENDENCY>_SOURCE_DIR cache variable
    set(_override_var "${_dep_upper}_SOURCE_DIR")
    if(DEFINED ${_override_var})
        get_filename_component(_candidate "${${_override_var}}" ABSOLUTE "${CMAKE_SOURCE_DIR}")
        if(IS_DIRECTORY "${_candidate}")
            set(_source_dir "${_candidate}")
            set(_use_local TRUE)
        else()
            message(FATAL_ERROR "${_override_var}='${_candidate}' is not a directory")
        endif()
    endif()

    if(_use_local)
        if("${_dep_lower}" STREQUAL "libneo")
            # SIMPLE controls its own tests with SIMPLE_TESTING. Keep libneo's
            # tests out of the parent target namespace, including older libneo
            # revisions that inherit BUILD_TESTING from the consumer.
            set(LIBNEO_BUILD_TESTING OFF CACHE BOOL "" FORCE)
            set(BUILD_TESTING OFF)
        endif()
        add_subdirectory("${_source_dir}" "${CMAKE_CURRENT_BINARY_DIR}/${_dep_lower}" EXCLUDE_FROM_ALL)
        return()
    endif()

    if("${_dep_lower}" STREQUAL "libneo")
        set(_repo_url "https://github.com/itpplasma/libneo.git")
        # LIBNEO_REF is the pin from cmake/dependency_versions.json unless
        # -DLIBNEO_REF=<branch|tag|sha> overrides it; never read from $ENV.
        string(STRIP "${LIBNEO_REF}" _branch)
        if(_branch STREQUAL "")
            message(FATAL_ERROR "LIBNEO_REF is empty; set it or pin libneo in cmake/dependency_versions.json")
        endif()

        # See the local-source case above. BUILD_TESTING is scoped to this
        # function so SIMPLE's own SIMPLE_TESTING setting is unaffected.
        set(LIBNEO_BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(BUILD_TESTING OFF)

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
