find_package(Python COMPONENTS Interpreter)

set(PYSIMPLE_PATH
    "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python")

# Register one pytest file as one CTest test named after the file.
function(add_pytest name labels timeout)
    add_test(NAME ${name}
        COMMAND ${Python_EXECUTABLE} -m pytest
            ${CMAKE_CURRENT_SOURCE_DIR}/${name}.py -v ${ARGN})
    set_tests_properties(${name} PROPERTIES
        ENVIRONMENT "${PYSIMPLE_PATH}"
        LABELS "python;${labels}"
        TIMEOUT ${timeout}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

add_pytest(test_simple_api "simple_api;slow" 900)
set_tests_properties(test_simple_api PROPERTIES
    ENVIRONMENT "${PYSIMPLE_PATH};SIMPLE_TEST_CHARTMAP=${CMAKE_BINARY_DIR}/test/tests/test_boozer_chartmap.nc"
    DEPENDS generate_test_boozer_chartmap_data)
add_pytest(test_orbit_macro "orbit;slow" 300)
add_pytest(test_classify_parity "classification;slow" 300)
add_pytest(test_examples "examples;slow" 300)
add_pytest(test_netcdf_results "netcdf;slow" 300)
add_pytest(test_chartmap_boundary_param "chartmap;slow" 300)
add_pytest(test_engineering "engineering" 300)
add_pytest(test_coil_clearance "engineering" 120)
add_pytest(test_plotting "plotting" 120)
add_pytest(test_cli "plotting" 120)
add_pytest(test_golden_orbit_results "golden" 120)
add_pytest(test_step_convergence "tools" 120)
