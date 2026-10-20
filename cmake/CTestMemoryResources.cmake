# One mem_gb slot represents one GiB of peak resident memory. Leave room for
# the OS, other jobs and tests without a memory reservation.
set(SIMPLE_TEST_MEMORY_GB "" CACHE STRING
    "Cap the CTest memory budget in GiB (empty: automatic)")
if(NOT SIMPLE_TEST_MEMORY_GB STREQUAL "" AND
   NOT SIMPLE_TEST_MEMORY_GB MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "SIMPLE_TEST_MEMORY_GB must be a positive integer or empty")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/CTestAvailableMemory.cmake")
simple_test_available_memory(_mem_bytes "/proc" "/sys/fs/cgroup")
math(EXPR _mem_gib "${_mem_bytes} / 1073741824")
math(EXPR _mem_reserve "(${_mem_gib} + 3) / 4")
if(_mem_reserve LESS 4)
    set(_mem_reserve 4)
endif()
math(EXPR _mem_slots "${_mem_gib} - ${_mem_reserve}")
if(_mem_slots LESS 0)
    set(_mem_slots 0)
endif()
if(NOT SIMPLE_TEST_MEMORY_GB STREQUAL "" AND
   SIMPLE_TEST_MEMORY_GB LESS _mem_slots)
    set(_mem_slots "${SIMPLE_TEST_MEMORY_GB}")
endif()

file(WRITE "${CMAKE_BINARY_DIR}/ctest-resources.json"
    "{\n  \"version\": {\"major\": 1, \"minor\": 0},\n"
    "  \"local\": [{\"mem_gb\": [{\"id\": \"ram\", \"slots\": ${_mem_slots}}]}]\n}\n")
message(STATUS "CTest memory budget: ${_mem_slots} GiB (${_mem_gib} GiB available)")

# CMake persists this in every generated CTestTestfile.cmake, so plain ctest
# and --test-dir both allocate memory. An explicit --resource-spec-file wins.
if(NOT CTEST_RESOURCE_SPEC_FILE)
    set(CTEST_RESOURCE_SPEC_FILE "${CMAKE_BINARY_DIR}/ctest-resources.json")
endif()
