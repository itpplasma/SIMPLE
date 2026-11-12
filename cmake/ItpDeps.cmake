# ItpDeps 1.0.0 (https://github.com/itpplasma/cmake-deps)
#
# Adds ITP projects and their shared dependencies to a CMake build, each once
# and pinned to a commit. The outermost project's pins and its copy of this
# file apply to the whole build. Copy this file unchanged from a release; see
# the cmake-deps README for the interface.
#
# SPDX-License-Identifier: MIT

set(_itp_deps_version 1.0.0)

# The first copy loaded serves the whole build. Later copies only compare
# versions, so a dependency cannot replace the outermost project's functions.
get_property(_itp_deps_active GLOBAL PROPERTY ITP_DEPS_VERSION)
if(_itp_deps_active)
    get_property(_itp_deps_active_file GLOBAL PROPERTY ITP_DEPS_FILE)
    string(REGEX MATCH "^[0-9]+" _itp_deps_new_major "${_itp_deps_version}")
    string(REGEX MATCH "^[0-9]+" _itp_deps_old_major "${_itp_deps_active}")
    if(_itp_deps_new_major GREATER _itp_deps_old_major)
        message(FATAL_ERROR
            "ItpDeps ${_itp_deps_version} at ${CMAKE_CURRENT_LIST_FILE} needs a "
            "newer major version than the active ItpDeps ${_itp_deps_active} at "
            "${_itp_deps_active_file}. Update the outermost project's copy.")
    elseif(_itp_deps_version VERSION_GREATER _itp_deps_active)
        set(_itp_deps_msg
            "ItpDeps ${_itp_deps_version} at ${CMAKE_CURRENT_LIST_FILE} is newer "
            "than the active ItpDeps ${_itp_deps_active} at ${_itp_deps_active_file}. "
            "Update the outermost project's copy.")
        if(ITP_DEPS_STRICT)
            message(FATAL_ERROR ${_itp_deps_msg})
        else()
            message(WARNING ${_itp_deps_msg})
        endif()
    endif()
    unset(_itp_deps_active)
    unset(_itp_deps_active_file)
    unset(_itp_deps_new_major)
    unset(_itp_deps_old_major)
    unset(_itp_deps_msg)
    unset(_itp_deps_version)
    return()
endif()
set_property(GLOBAL PROPERTY ITP_DEPS_VERSION "${_itp_deps_version}")
set_property(GLOBAL PROPERTY ITP_DEPS_FILE "${CMAKE_CURRENT_LIST_FILE}")
set_property(GLOBAL PROPERTY ITP_DEPS_NAMES "")
unset(_itp_deps_active)
unset(_itp_deps_version)

include(FetchContent)

function(_itp_deps_problem)
    list(JOIN ARGV "" message_text)
    if(ITP_DEPS_STRICT)
        message(FATAL_ERROR "${message_text}")
    else()
        message(WARNING "${message_text}")
    endif()
endfunction()

function(_itp_deps_git out_result out_output)
    set(_git git)
    if(GIT_EXECUTABLE)
        set(_git "${GIT_EXECUTABLE}")
    endif()
    execute_process(COMMAND "${_git}" -c advice.detachedHead=false ${ARGN}
        RESULT_VARIABLE _result
        OUTPUT_VARIABLE _output
        ERROR_VARIABLE _error
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(${out_result} "${_result}" PARENT_SCOPE)
    set(${out_output} "${_output}${_error}" PARENT_SCOPE)
endfunction()

function(_itp_deps_check_name name)
    if(NOT name MATCHES "^[a-z][a-z0-9_]*$")
        message(FATAL_ERROR
            "ItpDeps: dependency name '${name}' must match [a-z][a-z0-9_]*")
    endif()
endfunction()

function(_itp_deps_check_sha name sha what)
    if(NOT sha MATCHES "^[0-9a-f]+$")
        message(FATAL_ERROR
            "ItpDeps: ${what} for ${name} is '${sha}', not a full commit hash")
    endif()
    string(LENGTH "${sha}" _length)
    if(NOT _length EQUAL 40)
        message(FATAL_ERROR
            "ItpDeps: ${what} for ${name} is '${sha}', not a full 40-character "
            "commit hash")
    endif()
endfunction()

function(_itp_deps_url name explicit_url out_url)
    if(explicit_url)
        set(${out_url} "${explicit_url}" PARENT_SCOPE)
        return()
    endif()
    get_property(_pin_url GLOBAL PROPERTY ITP_DEPS_PIN_URL_${name})
    if(_pin_url)
        set(${out_url} "${_pin_url}" PARENT_SCOPE)
    else()
        set(${out_url} "https://github.com/itpplasma/${name}.git" PARENT_SCOPE)
    endif()
endfunction()

# Resolves a user override (commit, branch, tag or full ref) to a commit.
function(_itp_deps_resolve name url ref out_sha)
    if(ref MATCHES "^[0-9a-f]+$")
        string(LENGTH "${ref}" _length)
        if(_length EQUAL 40)
            set(${out_sha} "${ref}" PARENT_SCOPE)
            return()
        endif()
    endif()
    if($CACHE{FETCHCONTENT_FULLY_DISCONNECTED})
        message(FATAL_ERROR
            "ItpDeps: ITP_${name}_REF='${ref}' needs network access to resolve, "
            "but FETCHCONTENT_FULLY_DISCONNECTED is ON. Pass a full commit hash.")
    endif()
    _itp_deps_git(_result _output ls-remote "${url}" "${ref}")
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "ItpDeps: git ls-remote ${url} ${ref} failed:\n${_output}")
    endif()
    set(_sha "")
    string(REPLACE "\n" ";" _lines "${_output}")
    foreach(_candidate IN ITEMS "refs/heads/${ref}" "refs/tags/${ref}^{}"
            "refs/tags/${ref}" "${ref}")
        foreach(_line IN LISTS _lines)
            if(_line MATCHES "^([0-9a-f]+)\t(.*)$" AND CMAKE_MATCH_2 STREQUAL _candidate)
                set(_sha "${CMAKE_MATCH_1}")
                break()
            endif()
        endforeach()
        if(_sha)
            break()
        endif()
    endforeach()
    if(NOT _sha)
        message(FATAL_ERROR
            "ItpDeps: '${ref}' is neither a full commit hash nor a branch, tag or "
            "ref of ${url}")
    endif()
    set(${out_sha} "${_sha}" PARENT_SCOPE)
endfunction()

# Brings dir to commit sha without discarding anyone's work.
function(_itp_deps_checkout name url sha dir)
    if(NOT EXISTS "${dir}")
        if($CACHE{FETCHCONTENT_FULLY_DISCONNECTED})
            message(FATAL_ERROR
                "ItpDeps: ${dir} is missing and FETCHCONTENT_FULLY_DISCONNECTED is ON")
        endif()
        file(MAKE_DIRECTORY "${dir}")
        _itp_deps_git(_result _output -C "${dir}" init -q)
        if(NOT _result EQUAL 0)
            message(FATAL_ERROR "ItpDeps: git init in ${dir} failed:\n${_output}")
        endif()
    elseif(NOT EXISTS "${dir}/.git")
        message(FATAL_ERROR
            "ItpDeps: ${dir} exists but is not a git checkout. Remove it, or point "
            "FETCHCONTENT_SOURCE_DIR_${name} at the folder you want to use.")
    else()
        _itp_deps_git(_result _head -C "${dir}" rev-parse --verify -q HEAD)
        if(_result EQUAL 0 AND _head STREQUAL sha)
            return()
        endif()
        _itp_deps_git(_result _status -C "${dir}" status --porcelain
            --untracked-files=no)
        if(NOT _result EQUAL 0 OR NOT _status STREQUAL "")
            message(FATAL_ERROR
                "ItpDeps: ${dir} has local changes and is not at the requested "
                "commit ${sha}. Commit or remove them; ItpDeps never discards work.")
        endif()
    endif()
    _itp_deps_git(_result _output -C "${dir}" cat-file -e "${sha}^{commit}")
    if(NOT _result EQUAL 0)
        if($CACHE{FETCHCONTENT_FULLY_DISCONNECTED})
            message(FATAL_ERROR
                "ItpDeps: commit ${sha} of ${url} is not in ${dir} and "
                "FETCHCONTENT_FULLY_DISCONNECTED is ON")
        endif()
        message(STATUS "ItpDeps: fetching ${url} ${sha}")
        # Fetching a single commit by hash also reaches pull-request commits.
        _itp_deps_git(_result _output -C "${dir}" fetch -q --depth 1 "${url}" "${sha}")
        if(NOT _result EQUAL 0)
            message(FATAL_ERROR "ItpDeps: fetching ${url} ${sha} failed:\n${_output}")
        endif()
    endif()
    _itp_deps_git(_result _output -C "${dir}" checkout -q --detach "${sha}")
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "ItpDeps: checking out ${sha} in ${dir} failed:\n${_output}")
    endif()
endfunction()

function(_itp_deps_compare name caller requested)
    get_property(_how GLOBAL PROPERTY ITP_DEPS_HOW_${name})
    get_property(_used GLOBAL PROPERTY ITP_DEPS_COMMIT_${name})
    if(NOT _how MATCHES "^(pin|ref)$" OR requested STREQUAL _used)
        return()
    endif()
    get_property(_accept GLOBAL PROPERTY ITP_DEPS_PIN_ACCEPT_${name})
    if(requested IN_LIST _accept)
        return()
    endif()
    _itp_deps_problem(
        "ItpDeps: ${caller} requests ${name} ${requested}, but the build uses "
        "${_used} (chosen by ${_how}).")
endfunction()

function(itp_pin name sha)
    cmake_parse_arguments(PARSE_ARGV 2 _arg "" "GIT_REPOSITORY" "ACCEPT")
    if(_arg_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "ItpDeps: itp_pin(${name}): unknown arguments "
            "${_arg_UNPARSED_ARGUMENTS}")
    endif()
    _itp_deps_check_name("${name}")
    _itp_deps_check_sha("${name}" "${sha}" "pin")
    foreach(_accepted IN LISTS _arg_ACCEPT)
        _itp_deps_check_sha("${name}" "${_accepted}" "ACCEPT entry")
    endforeach()
    get_property(_pinned GLOBAL PROPERTY ITP_DEPS_PIN_${name})
    if(_pinned)
        if(NOT _pinned STREQUAL sha)
            get_property(_file GLOBAL PROPERTY ITP_DEPS_PIN_FILE_${name})
            _itp_deps_problem(
                "ItpDeps: ${CMAKE_CURRENT_LIST_FILE} pins ${name} ${sha}, but "
                "${_file} already pinned ${_pinned}, which wins.")
        endif()
        return()
    endif()
    get_property(_used GLOBAL PROPERTY ITP_DEPS_COMMIT_${name})
    if(_used)
        if(NOT _used STREQUAL sha)
            _itp_deps_problem(
                "ItpDeps: ${CMAKE_CURRENT_LIST_FILE} pins ${name} ${sha} after "
                "${name} was already added at ${_used}. Pin before adding.")
        endif()
        return()
    endif()
    set_property(GLOBAL PROPERTY ITP_DEPS_PIN_${name} "${sha}")
    set_property(GLOBAL PROPERTY ITP_DEPS_PIN_FILE_${name} "${CMAKE_CURRENT_LIST_FILE}")
    set_property(GLOBAL PROPERTY ITP_DEPS_PIN_ACCEPT_${name} "${_arg_ACCEPT}")
    if(_arg_GIT_REPOSITORY)
        set_property(GLOBAL PROPERTY ITP_DEPS_PIN_URL_${name} "${_arg_GIT_REPOSITORY}")
    endif()
endfunction()

function(itp_add_dependency name)
    cmake_parse_arguments(PARSE_ARGV 1 _arg "" "REF;GIT_REPOSITORY" "TARGET;OPTIONS")
    if(_arg_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "ItpDeps: itp_add_dependency(${name}): unknown "
            "arguments ${_arg_UNPARSED_ARGUMENTS}")
    endif()
    _itp_deps_check_name("${name}")
    _itp_deps_check_sha("${name}" "${_arg_REF}" "REF")
    string(TOUPPER "${name}" _upper)
    set(_caller "${CMAKE_CURRENT_LIST_FILE}")

    get_property(_callers GLOBAL PROPERTY ITP_DEPS_CALLERS_${name})
    set(_entry "${_caller}|${_arg_REF}")
    set(_repeat OFF)
    if(_entry IN_LIST _callers)
        set(_repeat ON)
    else()
        set_property(GLOBAL APPEND PROPERTY ITP_DEPS_CALLERS_${name} "${_entry}")
    endif()

    get_property(_used GLOBAL PROPERTY ITP_DEPS_COMMIT_${name})
    get_property(_how GLOBAL PROPERTY ITP_DEPS_HOW_${name})
    if(_how)
        if(NOT _repeat)
            _itp_deps_compare("${name}" "${_caller}" "${_arg_REF}")
        endif()
        return()
    endif()

    FetchContent_GetProperties(${name})
    set(_external OFF)
    if(${name}_POPULATED)
        set(_external ON)
    endif()
    foreach(_target IN LISTS _arg_TARGET)
        if(TARGET ${_target})
            set(_external ON)
        endif()
    endforeach()
    if(_external)
        set_property(GLOBAL APPEND PROPERTY ITP_DEPS_NAMES "${name}")
        set_property(GLOBAL PROPERTY ITP_DEPS_HOW_${name} external)
        set_property(GLOBAL PROPERTY ITP_DEPS_COMMIT_${name} "")
        set_property(GLOBAL PROPERTY ITP_DEPS_DIR_${name} "${${name}_SOURCE_DIR}")
        _itp_deps_problem(
            "ItpDeps: ${name} was added to the build without ItpDeps before "
            "${_caller} asked for it; its version is not checked.")
        return()
    endif()

    _itp_deps_url("${name}" "${_arg_GIT_REPOSITORY}" _url)
    get_property(_pin GLOBAL PROPERTY ITP_DEPS_PIN_${name})
    set(_dir "${FETCHCONTENT_BASE_DIR}/${name}-src")
    if(FETCHCONTENT_SOURCE_DIR_${_upper})
        set(_how local)
        set(_dir "${FETCHCONTENT_SOURCE_DIR_${_upper}}")
        _itp_deps_git(_result _sha -C "${_dir}" rev-parse --verify -q HEAD)
        if(NOT _result EQUAL 0)
            set(_sha "")
        endif()
    elseif(ITP_${_upper}_REF)
        set(_how override)
        _itp_deps_resolve("${_upper}" "${_url}" "${ITP_${_upper}_REF}" _sha)
        if(NOT _sha STREQUAL ITP_${_upper}_REF)
            set(_how unpinned)
            _itp_deps_problem(
                "ItpDeps: ${name} follows ITP_${_upper}_REF=${ITP_${_upper}_REF}, "
                "currently ${_sha}; the build is not reproducible.")
        endif()
    elseif(_pin)
        set(_how pin)
        set(_sha "${_pin}")
    else()
        set(_how ref)
        set(_sha "${_arg_REF}")
    endif()
    if(NOT _how STREQUAL "local")
        _itp_deps_checkout("${_upper}" "${_url}" "${_sha}" "${_dir}")
    endif()

    set_property(GLOBAL APPEND PROPERTY ITP_DEPS_NAMES "${name}")
    set_property(GLOBAL PROPERTY ITP_DEPS_HOW_${name} "${_how}")
    set_property(GLOBAL PROPERTY ITP_DEPS_COMMIT_${name} "${_sha}")
    set_property(GLOBAL PROPERTY ITP_DEPS_DIR_${name} "${_dir}")
    set_property(GLOBAL PROPERTY ITP_DEPS_URL_${name} "${_url}")
    _itp_deps_compare("${name}" "${_caller}" "${_arg_REF}")

    # Settings below stay in this function: the dependency sees them, the
    # caller and the cache do not.
    set(BUILD_TESTING OFF)
    foreach(_option IN LISTS _arg_OPTIONS)
        if(NOT _option MATCHES "^([A-Za-z_][A-Za-z0-9_]*) (.*)$")
            message(FATAL_ERROR
                "ItpDeps: OPTIONS entry '${_option}' for ${name} is not '<VAR> <value>'")
        endif()
        set(${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
    endforeach()
    # The checkout is ready, so FetchContent only adds it. FetchContent copies
    # this variable into the cache, where it would look like a user override
    # on the next configure; remove that entry again.
    if(DEFINED CACHE{FETCHCONTENT_SOURCE_DIR_${_upper}})
        set(_user_entry ON)
    else()
        set(_user_entry OFF)
    endif()
    set(FETCHCONTENT_SOURCE_DIR_${_upper} "${_dir}")
    FetchContent_Declare(${name} GIT_REPOSITORY "${_url}" GIT_TAG "${_sha}")
    FetchContent_MakeAvailable(${name})
    if(NOT _user_entry)
        unset(FETCHCONTENT_SOURCE_DIR_${_upper} CACHE)
    endif()
    if(EXISTS "${_dir}/CMakeLists.txt")
        # Build only what the caller links; set after adding so this also
        # works before CMake 3.28, whose FetchContent_Declare lacks the option.
        set_property(DIRECTORY "${_dir}" PROPERTY EXCLUDE_FROM_ALL TRUE)
    endif()
endfunction()

function(_itp_deps_json_string value out_var)
    string(REPLACE "\\" "\\\\" value "${value}")
    string(REPLACE "\"" "\\\"" value "${value}")
    set(${out_var} "\"${value}\"" PARENT_SCOPE)
endfunction()

function(_itp_deps_report)
    get_property(_names GLOBAL PROPERTY ITP_DEPS_NAMES)
    get_property(_version GLOBAL PROPERTY ITP_DEPS_VERSION)
    get_property(_file GLOBAL PROPERTY ITP_DEPS_FILE)
    if(NOT _names)
        return()
    endif()
    message(STATUS "ItpDeps ${_version} (${_file}):")
    _itp_deps_json_string("${_version}" _jversion)
    _itp_deps_json_string("${_file}" _jfile)
    set(_json "{\n  \"itp_deps_version\": ${_jversion},\n  \"itp_deps_file\": ${_jfile},\n  \"dependencies\": [")
    set(_separator "")
    foreach(_name IN LISTS _names)
        get_property(_sha GLOBAL PROPERTY ITP_DEPS_COMMIT_${_name})
        get_property(_how GLOBAL PROPERTY ITP_DEPS_HOW_${_name})
        get_property(_dir GLOBAL PROPERTY ITP_DEPS_DIR_${_name})
        get_property(_url GLOBAL PROPERTY ITP_DEPS_URL_${_name})
        get_property(_callers GLOBAL PROPERTY ITP_DEPS_CALLERS_${_name})
        set(_shown "${_sha}")
        if(NOT _shown)
            set(_shown "unknown")
        endif()
        message(STATUS "  ${_name} ${_shown} (${_how})")
        set(_jcallers "")
        set(_caller_separator "")
        foreach(_entry IN LISTS _callers)
            string(REGEX MATCH "^(.*)\\|([0-9a-f]+)$" _match "${_entry}")
            _itp_deps_json_string("${CMAKE_MATCH_1}" _jcaller)
            _itp_deps_json_string("${CMAKE_MATCH_2}" _jref)
            string(APPEND _jcallers
                "${_caller_separator}\n        {\"file\": ${_jcaller}, \"ref\": ${_jref}}")
            set(_caller_separator ",")
        endforeach()
        foreach(_field IN ITEMS name sha how dir url)
            if(_field STREQUAL "name")
                _itp_deps_json_string("${_name}" _j${_field})
            else()
                _itp_deps_json_string("${_${_field}}" _j${_field})
            endif()
        endforeach()
        string(APPEND _json "${_separator}\n    {\"name\": ${_jname}, \"commit\": ${_jsha}, "
            "\"how\": ${_jhow},\n      \"source_dir\": ${_jdir}, \"url\": ${_jurl},\n"
            "      \"callers\": [${_jcallers}]}")
        set(_separator ",")
    endforeach()
    string(APPEND _json "\n  ]\n}\n")
    file(WRITE "${CMAKE_BINARY_DIR}/itp-deps.json" "${_json}")
endfunction()

cmake_language(DEFER DIRECTORY "${CMAKE_SOURCE_DIR}" CALL _itp_deps_report)
