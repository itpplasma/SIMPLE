add_test(NAME test_ctest_memory_resources COMMAND ${CMAKE_COMMAND}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ctest_memory_resources.cmake)
set_tests_properties(test_ctest_memory_resources PROPERTIES LABELS "unit" TIMEOUT 30)

add_subdirectory(tests)
if(ENABLE_PYTHON_INTERFACE AND SIMPLE_PYTHON_INTERFACE_BUILT)
    add_subdirectory(python)
elseif(ENABLE_PYTHON_INTERFACE)
    message(STATUS "Skipping Python tests because the SIMPLE Python backend was not built")
endif()
