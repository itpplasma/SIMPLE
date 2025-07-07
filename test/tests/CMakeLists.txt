##################################################
# DOWNLOAD TEST DATA
##################################################
file(DOWNLOAD
    https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc
    ${CMAKE_CURRENT_BINARY_DIR}/wout.nc
)

# Convert VMEC to GVEC format for testing field_vmec vs field_gvec
# This requires GVEC Python API or command-line tools
find_program(PYTHON_EXECUTABLE python3)
if(PYTHON_EXECUTABLE)
    # Try to convert VMEC to GVEC using Python
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "
try:
    import gvec
    print('GVEC Python module found, converting VMEC to GVEC format...')
    gvec.read_vmec('${CMAKE_CURRENT_BINARY_DIR}/wout.nc')
    gvec.save('${CMAKE_CURRENT_BINARY_DIR}/wout_gvec.dat')
    print('Successfully converted VMEC to GVEC format: wout_gvec.dat')
except ImportError:
    print('GVEC Python module not found, skipping VMEC to GVEC conversion')
except Exception as e:
    print(f'Failed to convert VMEC to GVEC: {e}')
"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        RESULT_VARIABLE GVEC_CONVERT_RESULT
    )
else()
    message(STATUS "Python not found, skipping VMEC to GVEC conversion")
endif()

##################################################
# INTERNAL UNIT TEST STUFF
##################################################
add_library(testing_utils ${CMAKE_CURRENT_SOURCE_DIR}/utility.f90 ${CMAKE_CURRENT_SOURCE_DIR}/results.f90)

##################################################
# TEST EXECUTABLES
##################################################
add_executable(utility_test.x utility_test.f90)
target_link_libraries(utility_test.x testing_utils)

if (ENABLE_OPENMP)
    add_executable (test_sympl_tok.x test_sympl.f90)
    target_link_libraries(test_sympl_tok.x simple)
    add_executable (test_sympl.x test_sympl_stell.f90)
    target_link_libraries(test_sympl.x simple)
endif()

add_executable (test_coord_trans.x test_coord_trans.f90)
target_link_libraries(test_coord_trans.x simple)

add_executable (test_magfie.x test_magfie.f90)
target_link_libraries(test_magfie.x simple)
add_executable (test_poincare1.x test_poincare1.f90)
target_link_libraries(test_poincare1.x simple)
add_executable (test_poiplot_classification.x test_poiplot_classification.f90
)
target_link_libraries(test_poiplot_classification.x simple)
add_executable (test_parse_ants.x test_parse_ants.f90)
target_link_libraries(test_parse_ants.x simple)
add_executable (test_boozer.x test_boozer.f90)
target_link_libraries(test_boozer.x simple)
add_executable (test_orbits.x test_orbits_vmec.f90)
target_link_libraries(test_orbits.x simple)
add_executable (test_collis.x ../../src/test_collis.f90)
target_link_libraries(test_collis.x simple)

add_executable (test_samplers.x test_samplers.f90)
target_link_libraries(test_samplers.x simple)

add_executable (test_coordinates.x test_coordinates.f90)
target_link_libraries(test_coordinates.x simple)
add_test(NAME test_coordinates COMMAND test_coordinates.x vmec cyl)

# Generate GVEC test data file for test_gvec using elliptic tokamak example
set(GVEC_TEST_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/../../build/_deps/gvec-src/test-CI/examples/analytic_gs_elliptok/parameter.ini")
set(GVEC_TEST_STATE "${CMAKE_CURRENT_SOURCE_DIR}/../../test/test_data/GVEC_elliptok_State_final.dat")

# Disable GVEC test data generation for now - causes build issues
# TODO: Fix GVEC test data generation path
# add_custom_command(
#     OUTPUT ${GVEC_TEST_STATE}
#     COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/../../test/test_data
#     COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_BINARY_DIR}
#             ${CMAKE_CURRENT_BINARY_DIR}/../../_deps/gvec-build/bin/gvec ${GVEC_TEST_INPUT}
#     COMMAND ${CMAKE_COMMAND} -E copy_if_different 
#             ${CMAKE_CURRENT_BINARY_DIR}/GVEC_SOLOVIEV_ELLIPTOK_State_0000_00000000.dat 
#             ${GVEC_TEST_STATE}
#     DEPENDS gvec
#     COMMENT "Generating GVEC elliptic tokamak test state file"
#     VERBATIM
# )

# Custom target to ensure the state file is generated
# add_custom_target(gvec_test_data DEPENDS ${GVEC_TEST_STATE})

add_executable (test_gvec.x test_gvec.f90)
target_link_libraries(test_gvec.x simple)
# add_dependencies(test_gvec.x gvec_test_data)
add_test(NAME test_gvec COMMAND test_gvec.x)

add_executable (test_vmec_gvec.x test_vmec_gvec.f90)
target_link_libraries(test_vmec_gvec.x simple)
add_test(NAME test_vmec_gvec COMMAND test_vmec_gvec.x)
set_tests_properties(test_vmec_gvec PROPERTIES
    LABELS "slow"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable (test_canonical_gvec.x test_canonical_gvec.f90)
target_link_libraries(test_canonical_gvec.x simple)
add_test(NAME test_canonical_gvec COMMAND test_canonical_gvec.x)

add_executable (test_simple_vmec_gvec.x test_simple_vmec_gvec.f90)
target_link_libraries(test_simple_vmec_gvec.x simple)
add_test(NAME test_simple_vmec_gvec COMMAND test_simple_vmec_gvec.x)
set_tests_properties(test_simple_vmec_gvec PROPERTIES
    LABELS "slow;integration"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable (test_adapter_consistency.x test_adapter_consistency.f90)
target_link_libraries(test_adapter_consistency.x simple)
add_test(NAME test_adapter_consistency COMMAND test_adapter_consistency.x)
set_tests_properties(test_adapter_consistency PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable (test_vmec_gvec_adapter.x test_vmec_gvec_adapter.f90)
target_link_libraries(test_vmec_gvec_adapter.x simple)
add_test(NAME test_vmec_gvec_adapter COMMAND test_vmec_gvec_adapter.x)
set_tests_properties(test_vmec_gvec_adapter PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LABELS "adapter"
)

add_executable (test_spline_3d.x test_spline_3d.f90)
target_link_libraries(test_spline_3d.x simple)
add_test(NAME test_spline_3d COMMAND test_spline_3d.x)

# Add 2D field export utility (not a test, but a diagnostic tool)
add_executable (export_field_2d.x export_field_2d.f90)
target_link_libraries(export_field_2d.x simple)

add_test(NAME test_magfie COMMAND test_magfie.x)

# Unit tests for refactored modules
add_executable(test_stencil_utils.x test_stencil_utils.f90)
target_link_libraries(test_stencil_utils.x simple)
add_test(NAME test_stencil_utils COMMAND test_stencil_utils.x)

add_executable(test_array_utils.x test_array_utils.f90)
target_link_libraries(test_array_utils.x simple)
add_test(NAME test_array_utils COMMAND test_array_utils.x)

# Regression tests for canonical coordinates
# Note: The following test files were removed as they don't exist yet:
# - test_canonical_stencil.f90
# - test_canonical_full_regression.f90  
# - test_canonical_baseline_simple.f90

add_subdirectory(field_can)
add_subdirectory(magfie)
add_subdirectory(orbit)

##################################################
# GOLDEN RECORD TESTS
##################################################

# Copy test cases to build directory for golden record tests
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record/boozer 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/golden_record/test_cases/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record/canonical 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/golden_record/test_cases/)

# Golden record test: expected success (46f1f53 vs itself should match)
add_test(
    NAME golden_record_success_46f1f53
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record/golden_record.sh 46f1f53
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(golden_record_success_46f1f53 PROPERTIES
    TIMEOUT 300  # 5 minutes
    LABELS "golden_record;regression;slow"
    ENVIRONMENT "GOLDEN_RECORD_BASE_DIR=${CMAKE_CURRENT_BINARY_DIR}/golden_record"
)

# Golden record test: expected failure (current vs old version should not match)
add_test(
    NAME golden_record_failure
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record/golden_record.sh e25ca7e
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(golden_record_failure PROPERTIES
    TIMEOUT 300  # 5 minutes
    WILL_FAIL TRUE  # This test is expected to fail
    LABELS "golden_record;regression;slow"
    ENVIRONMENT "GOLDEN_RECORD_BASE_DIR=${CMAKE_CURRENT_BINARY_DIR}/golden_record"
)
