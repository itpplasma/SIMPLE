##################################################
# DOWNLOAD TEST DATA
##################################################
# Download VMEC test file to unified test_data directory (shared with Python tests)
set(TEST_DATA_DIR "${CMAKE_SOURCE_DIR}/test/test_data")
set(WOUT_FILE "${TEST_DATA_DIR}/wout.nc")
set(WOUT_URL "https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc")

# Download if not present
if(NOT EXISTS "${WOUT_FILE}")
    file(MAKE_DIRECTORY "${TEST_DATA_DIR}")
    message(STATUS "Downloading VMEC test file to ${WOUT_FILE}...")
    file(DOWNLOAD "${WOUT_URL}" "${WOUT_FILE}")
    message(STATUS "Downloaded VMEC test file")
endif()

# Create symlink in test binary directory for backward compatibility
file(CREATE_LINK "${WOUT_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/wout.nc" SYMBOLIC)

# Generate chartmap file from VMEC for testing chartmap coordinate detection
# Uses libneo's vmec_to_chartmap tool built as part of dependencies
set(CHARTMAP_FILE "${CMAKE_CURRENT_BINARY_DIR}/wout.chartmap.nc")

add_custom_command(
    OUTPUT ${CHARTMAP_FILE}
    COMMAND $<TARGET_FILE:vmec_to_chartmap.x> ${WOUT_FILE} ${CHARTMAP_FILE} 17 33 33
    DEPENDS ${WOUT_FILE} vmec_to_chartmap.x
    COMMENT "Generating chartmap test file from VMEC"
    VERBATIM
)
add_custom_target(chartmap_test_data DEPENDS ${CHARTMAP_FILE})

# map2disc-based chartmap generation (via libneo Python API)
if(SIMPLE_ENABLE_PYTHON_TOOLS)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    set(CHARTMAP_MAP2DISC_FILE "${CMAKE_CURRENT_BINARY_DIR}/wout.chartmap.map2disc.nc")
    set(CHARTMAP_NCSX_MAP2DISC_FILE "${CMAKE_CURRENT_BINARY_DIR}/wout_ncsx.chartmap.map2disc.nc")
    set(CHARTMAP_NCSX_MAP2DISC_WALL_FILE "${CMAKE_CURRENT_BINARY_DIR}/wout_ncsx.chartmap.map2disc.wall0p1m.nc")

    set(LIBNEO_PYTHON_DIR "" CACHE PATH "Path to libneo/python directory containing the libneo Python package")
    if(LIBNEO_PYTHON_DIR STREQUAL "")
        if(DEFINED libneo_SOURCE_DIR AND EXISTS "${libneo_SOURCE_DIR}/python/libneo/chartmap.py")
            set(LIBNEO_PYTHON_DIR "${libneo_SOURCE_DIR}/python")
        elseif(DEFINED LIBNEO_SOURCE_DIR AND EXISTS "${LIBNEO_SOURCE_DIR}/python/libneo/chartmap.py")
            set(LIBNEO_PYTHON_DIR "${LIBNEO_SOURCE_DIR}/python")
        endif()

        if(NOT LIBNEO_PYTHON_DIR STREQUAL "")
            set(LIBNEO_PYTHON_DIR "${LIBNEO_PYTHON_DIR}" CACHE PATH
                "Path to libneo/python directory containing the libneo Python package" FORCE)
        endif()
    endif()

    if(LIBNEO_PYTHON_DIR STREQUAL "")
        message(FATAL_ERROR
            "LIBNEO_PYTHON_DIR is empty but map2disc-based chartmap tests are enabled. "
            "Set -DLIBNEO_PYTHON_DIR=/path/to/libneo/python or disable via -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF.")
    endif()

    if(NOT EXISTS "${LIBNEO_PYTHON_DIR}/libneo/chartmap.py")
        message(FATAL_ERROR
            "libneo Python chartmap API not found at '${LIBNEO_PYTHON_DIR}/libneo/chartmap.py'. "
            "Set -DLIBNEO_PYTHON_DIR=/path/to/libneo/python or disable via -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF.")
    endif()

    add_test(NAME generate_chartmap_map2disc
        COMMAND ${CMAKE_COMMAND} -E env
            PYTHONPATH=${LIBNEO_PYTHON_DIR}
            ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/generate_chartmap_map2disc.py
            --libneo-python ${LIBNEO_PYTHON_DIR}
            --wout ${CMAKE_CURRENT_BINARY_DIR}/wout.nc
            --out ${CHARTMAP_MAP2DISC_FILE}
            --nrho 17 --ntheta 33 --nzeta 33
            --s-boundary 1.0
            --boundary-offset 0.0
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(generate_chartmap_map2disc PROPERTIES
        LABELS "integration;python"
        TIMEOUT 600
        FIXTURES_SETUP chartmap_map2disc)

    add_test(NAME generate_chartmap_ncsx_map2disc
        COMMAND ${CMAKE_COMMAND} -E env
            PYTHONPATH=${LIBNEO_PYTHON_DIR}
            ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/generate_chartmap_map2disc.py
            --libneo-python ${LIBNEO_PYTHON_DIR}
            --wout ${CMAKE_CURRENT_BINARY_DIR}/wout_ncsx.nc
            --out ${CHARTMAP_NCSX_MAP2DISC_FILE}
            --nrho 17 --ntheta 33 --nzeta 33
            --s-boundary 1.0
            --boundary-offset 0.0
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(generate_chartmap_ncsx_map2disc PROPERTIES
        LABELS "integration;python"
        TIMEOUT 600
        FIXTURES_SETUP chartmap_ncsx_map2disc)

    add_test(NAME generate_chartmap_ncsx_map2disc_wall
        COMMAND ${CMAKE_COMMAND} -E env
            PYTHONPATH=${LIBNEO_PYTHON_DIR}
            ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/generate_chartmap_map2disc.py
            --libneo-python ${LIBNEO_PYTHON_DIR}
            --wout ${CMAKE_CURRENT_BINARY_DIR}/wout_ncsx.nc
            --out ${CHARTMAP_NCSX_MAP2DISC_WALL_FILE}
            --nrho 17 --ntheta 33 --nzeta 33
            --s-boundary 1.0
            --boundary-offset 0.1
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(generate_chartmap_ncsx_map2disc_wall PROPERTIES
        LABELS "integration;python"
        TIMEOUT 600
        FIXTURES_SETUP chartmap_ncsx_map2disc_wall)

    add_test(NAME test_chartmap_ncsx_map2disc_wall_offset
        COMMAND ${CMAKE_COMMAND} -E env
            PYTHONPATH=${LIBNEO_PYTHON_DIR}
            ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/plot_chartmap_wall_offset_check.py
            --wout ${CMAKE_CURRENT_BINARY_DIR}/wout_ncsx.nc
            --chartmap ${CHARTMAP_NCSX_MAP2DISC_WALL_FILE}
            --out ${CMAKE_CURRENT_BINARY_DIR}/chartmap_wall_offset_check.png
            --offset 0.1
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(test_chartmap_ncsx_map2disc_wall_offset PROPERTIES
        LABELS "integration;plot;python"
        TIMEOUT 180
        FIXTURES_REQUIRED chartmap_ncsx_map2disc_wall)

    add_test(NAME plot_chartmap_vmec_vs_map2disc
        COMMAND ${CMAKE_COMMAND} -E env
            PYTHONPATH=${LIBNEO_PYTHON_DIR}
            ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/plot_chartmap_compare_map2disc.py
            --vmec ${CHARTMAP_FILE}
            --map2disc ${CHARTMAP_MAP2DISC_FILE}
            --out ${CMAKE_CURRENT_BINARY_DIR}/chartmap_vmec_vs_map2disc_slices.png
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(plot_chartmap_vmec_vs_map2disc PROPERTIES
        DEPENDS "chartmap_test_data;generate_chartmap_map2disc"
        LABELS "integration;plot;python"
        TIMEOUT 120
        FIXTURES_REQUIRED chartmap_map2disc)
endif()

# Download NCSX test data for field equivalence tests (coils + wout)
set(WOUT_NCSX_FILE "${TEST_DATA_DIR}/wout_ncsx.nc")
set(WOUT_NCSX_URL "https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc")
set(COILS_STELLOPT_FILE "${TEST_DATA_DIR}/coils.c09r00")
set(COILS_STELLOPT_URL "https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00")

if(NOT EXISTS "${WOUT_NCSX_FILE}")
    message(STATUS "Downloading NCSX VMEC file to ${WOUT_NCSX_FILE}...")
    file(DOWNLOAD "${WOUT_NCSX_URL}" "${WOUT_NCSX_FILE}")
endif()

if(NOT EXISTS "${COILS_STELLOPT_FILE}")
    message(STATUS "Downloading NCSX coils file to ${COILS_STELLOPT_FILE}...")
    file(DOWNLOAD "${COILS_STELLOPT_URL}" "${COILS_STELLOPT_FILE}")
endif()

# Generate chartmap from NCSX VMEC and convert coils for field equivalence test
set(CHARTMAP_NCSX_FILE "${CMAKE_CURRENT_BINARY_DIR}/wout_ncsx.chartmap.nc")
set(COILS_SIMPLE_FILE "${CMAKE_CURRENT_BINARY_DIR}/coils.simple")

add_custom_command(
    OUTPUT ${CHARTMAP_NCSX_FILE}
    COMMAND $<TARGET_FILE:vmec_to_chartmap.x> ${WOUT_NCSX_FILE} ${CHARTMAP_NCSX_FILE} 17 33 33
    DEPENDS ${WOUT_NCSX_FILE} vmec_to_chartmap.x
    COMMENT "Generating NCSX chartmap test file from VMEC"
    VERBATIM
)

# Convert STELLOPT coils format to SIMPLE format using Python
if(SIMPLE_ENABLE_PYTHON_TOOLS)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_custom_command(
        OUTPUT ${COILS_SIMPLE_FILE}
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/convert_coils_stellopt_to_simple.py
            ${COILS_STELLOPT_FILE} ${COILS_SIMPLE_FILE}
        DEPENDS ${COILS_STELLOPT_FILE}
            ${CMAKE_CURRENT_SOURCE_DIR}/convert_coils_stellopt_to_simple.py
        COMMENT "Converting STELLOPT coils to SIMPLE format"
        VERBATIM
    )
else()
    add_executable(convert_coils_stellopt_to_simple.x
        ${CMAKE_CURRENT_SOURCE_DIR}/convert_coils_stellopt_to_simple.f90)
    add_custom_command(
        OUTPUT ${COILS_SIMPLE_FILE}
        COMMAND $<TARGET_FILE:convert_coils_stellopt_to_simple.x>
            ${COILS_STELLOPT_FILE} ${COILS_SIMPLE_FILE}
        DEPENDS convert_coils_stellopt_to_simple.x ${COILS_STELLOPT_FILE}
        COMMENT "Converting STELLOPT coils to SIMPLE format (Fortran)"
        VERBATIM
    )
endif()

add_custom_target(ncsx_test_data DEPENDS ${CHARTMAP_NCSX_FILE} ${COILS_SIMPLE_FILE})

# Symlinks for NCSX test data
file(CREATE_LINK "${WOUT_NCSX_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/wout_ncsx.nc" SYMBOLIC)

# Convert VMEC to GVEC format for testing field_vmec vs field_gvec
# This requires GVEC Python API or command-line tools
if(ENABLE_GVEC)
    find_program(PYTHON_EXECUTABLE python3)
    if(PYTHON_EXECUTABLE AND SIMPLE_ENABLE_PYTHON_TOOLS)
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
endif()

##################################################
# INTERNAL UNIT TEST STUFF
##################################################
add_library(testing_utils ${CMAKE_CURRENT_SOURCE_DIR}/utility.f90 ${CMAKE_CURRENT_SOURCE_DIR}/results.f90)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_compile_options(testing_utils PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>)
else()
    target_compile_options(testing_utils PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)
endif()

# Link coverage flags to testing utilities when coverage is enabled
if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE MATCHES "Debug|Profile")
    target_link_libraries(testing_utils PUBLIC coverage_flags)
endif()

##################################################
# TEST EXECUTABLES
##################################################
add_executable(utility_test.x utility_test.f90)
target_link_libraries(utility_test.x testing_utils)

if (ENABLE_OPENMP)
    add_executable (test_sympl_tok.x test_sympl.f90)
    target_link_libraries(test_sympl_tok.x simple)
    add_test(NAME test_sympl_tokamak
        COMMAND test_sympl_tok.x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(test_sympl_tokamak PROPERTIES
        LABELS "integration")

    if(SIMPLE_ENABLE_PYTHON_TOOLS)
        add_test(NAME test_sympl_tokamak_plot
            COMMAND ${CMAKE_COMMAND} -E env python3
                ${CMAKE_CURRENT_SOURCE_DIR}/plot_tokamak_testfield.py
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(test_sympl_tokamak_plot PROPERTIES
            DEPENDS test_sympl_tokamak
            LABELS "integration;plot;python")
    endif()

    add_executable (test_sympl.x test_sympl_stell.f90)
    target_link_libraries(test_sympl.x simple)
    add_executable (test_sympl_testfield.x test_sympl_testfield.f90)
    target_link_libraries(test_sympl_testfield.x simple)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
        target_compile_options(test_sympl_testfield.x PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>)
    else()
        target_compile_options(test_sympl_testfield.x PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)
    endif()
    target_compile_definitions(test_sympl_testfield.x PRIVATE WOUT_FILE="${WOUT_FILE}")
    add_test(NAME test_sympl_testfield
        COMMAND test_sympl_testfield.x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
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

# GVEC-related tests (only if GVEC support is enabled)
if(ENABLE_GVEC)
    add_executable (test_gvec.x test_gvec.f90)
    target_link_libraries(test_gvec.x simple)
    # add_dependencies(test_gvec.x gvec_test_data)
    add_test(NAME test_gvec COMMAND test_gvec.x)

    add_executable (test_vmec_gvec.x test_vmec_gvec.f90)
    target_link_libraries(test_vmec_gvec.x simple)
    add_test(NAME test_vmec_gvec COMMAND test_vmec_gvec.x)

    add_executable (test_canonical_gvec.x test_canonical_gvec.f90)
    target_link_libraries(test_canonical_gvec.x simple)
    add_test(NAME test_canonical_gvec COMMAND test_canonical_gvec.x)

    add_executable (test_simple_vmec_gvec.x test_simple_vmec_gvec.f90)
    target_link_libraries(test_simple_vmec_gvec.x simple)
    add_test(NAME test_simple_vmec_gvec COMMAND test_simple_vmec_gvec.x)
    set_tests_properties(test_simple_vmec_gvec PROPERTIES
        LABELS "slow;integration"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        WILL_FAIL TRUE  # Expected to fail
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
endif()


# Add 2D field export utility (requires GVEC support)
if(ENABLE_GVEC)
    add_executable (export_field_2d.x export_field_2d.f90)
    target_link_libraries(export_field_2d.x simple)
endif()

add_test(NAME test_magfie COMMAND test_magfie.x)

# Unit tests for refactored modules
add_executable(test_stencil_utils.x test_stencil_utils.f90)
target_link_libraries(test_stencil_utils.x simple)
add_test(NAME test_stencil_utils COMMAND test_stencil_utils.x)

add_executable(test_array_utils.x test_array_utils.f90)
target_link_libraries(test_array_utils.x simple)
add_test(NAME test_array_utils COMMAND test_array_utils.x)

# New comprehensive unit tests
add_executable(test_util_simple.x test_util_simple.f90)
target_link_libraries(test_util_simple.x simple)
add_test(NAME test_util_simple COMMAND test_util_simple.x)

add_executable(test_timing.x test_timing.f90)
target_link_libraries(test_timing.x simple)
add_test(NAME test_timing COMMAND test_timing.x)

add_executable(test_sorting.x test_sorting.f90)
target_link_libraries(test_sorting.x simple)
add_test(NAME test_sorting COMMAND test_sorting.x)

add_executable(test_lapack_interfaces.x test_lapack_interfaces.f90)
target_link_libraries(test_lapack_interfaces.x simple)
add_test(NAME test_lapack_interfaces COMMAND test_lapack_interfaces.x)

add_executable(test_orbit_symplectic_base.x test_orbit_symplectic_base.f90)
target_link_libraries(test_orbit_symplectic_base.x simple)
add_test(NAME test_orbit_symplectic_base COMMAND test_orbit_symplectic_base.x)

add_executable(test_field_base.x test_field_base.f90)
target_link_libraries(test_field_base.x simple)
add_test(NAME test_field_base COMMAND test_field_base.x)

add_executable(test_coordinates_simple.x test_coordinates_simple.f90)
target_link_libraries(test_coordinates_simple.x simple)
add_test(NAME test_coordinates_simple COMMAND test_coordinates_simple.x)

add_executable(test_coordinate_scaling.x test_coordinate_scaling.f90)
target_link_libraries(test_coordinate_scaling.x simple)
add_test(NAME test_coordinate_scaling COMMAND test_coordinate_scaling.x)
set_tests_properties(test_coordinate_scaling PROPERTIES LABELS "unit")

add_executable(test_coordinate_refactoring.x test_coordinate_refactoring.f90)
target_link_libraries(test_coordinate_refactoring.x simple)
add_test(NAME test_coordinate_refactoring
    COMMAND test_coordinate_refactoring.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_coordinate_refactoring PROPERTIES
    LABELS "integration")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set_tests_properties(test_coordinate_refactoring PROPERTIES WILL_FAIL TRUE)
endif()

add_executable(test_splined_field_derivatives.x test_splined_field_derivatives.f90)
target_link_libraries(test_splined_field_derivatives.x simple)
add_dependencies(test_splined_field_derivatives.x ncsx_test_data)
add_test(NAME test_splined_field_derivatives
    COMMAND test_splined_field_derivatives.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_splined_field_derivatives PROPERTIES
    LABELS "integration")

add_executable(test_bder_from_splines.x test_bder_from_splines.f90)
target_link_libraries(test_bder_from_splines.x simple)
add_dependencies(test_bder_from_splines.x ncsx_test_data)
add_test(NAME test_bder_from_splines
    COMMAND test_bder_from_splines.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_bder_from_splines PROPERTIES
    LABELS "integration")

add_executable(test_refcoords_file_detection.x test_refcoords_file_detection.f90)
target_link_libraries(test_refcoords_file_detection.x simple)
add_dependencies(test_refcoords_file_detection.x chartmap_test_data)
add_test(NAME test_refcoords_file_detection
    COMMAND test_refcoords_file_detection.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_refcoords_file_detection PROPERTIES
    LABELS "unit")

add_executable(test_field_equivalence_chartmap.x test_field_equivalence_chartmap.f90)
target_link_libraries(test_field_equivalence_chartmap.x simple)
add_dependencies(test_field_equivalence_chartmap.x ncsx_test_data)

if(SIMPLE_ENABLE_PYTHON_TOOLS)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/plot_field_equivalence.py
        ${CMAKE_CURRENT_BINARY_DIR}/plot_field_equivalence.py
        COPYONLY)
endif()

add_test(NAME test_field_equivalence_chartmap
    COMMAND test_field_equivalence_chartmap.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_field_equivalence_chartmap PROPERTIES
    LABELS "integration"
    TIMEOUT 120)

if(SIMPLE_ENABLE_PYTHON_TOOLS)
    add_test(NAME test_field_equivalence_chartmap_map2disc
        COMMAND test_field_equivalence_chartmap.x wout_ncsx.chartmap.map2disc.nc field_equiv_map2disc
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(test_field_equivalence_chartmap_map2disc PROPERTIES
        LABELS "integration;python"
        TIMEOUT 120
        DEPENDS "generate_chartmap_ncsx_map2disc"
        FIXTURES_REQUIRED chartmap_ncsx_map2disc)

    add_test(NAME test_field_equivalence_chartmap_map2disc_wall
        COMMAND test_field_equivalence_chartmap.x wout_ncsx.chartmap.map2disc.wall0p1m.nc field_equiv_map2disc_wall
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(test_field_equivalence_chartmap_map2disc_wall PROPERTIES
        LABELS "integration;python"
        TIMEOUT 120
        DEPENDS "generate_chartmap_ncsx_map2disc_wall"
        FIXTURES_REQUIRED chartmap_ncsx_map2disc_wall)
endif()

# Chartmap pipeline unit tests (issue #227)
add_executable(test_chartmap_pipeline.x test_chartmap_pipeline.f90)
target_link_libraries(test_chartmap_pipeline.x simple)
add_dependencies(test_chartmap_pipeline.x ncsx_test_data)

if(SIMPLE_ENABLE_PYTHON_TOOLS)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/plot_chartmap_pipeline.py
        ${CMAKE_CURRENT_BINARY_DIR}/plot_chartmap_pipeline.py
        COPYONLY)
endif()

add_test(NAME test_chartmap_pipeline
    COMMAND test_chartmap_pipeline.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_chartmap_pipeline PROPERTIES
    LABELS "unit"
    TIMEOUT 60)

# Batch spline optimization test
add_executable(test_batch_splines.x ../test_batch_splines.f90)
target_link_libraries(test_batch_splines.x simple testing_utils)
add_test(NAME test_batch_splines COMMAND test_batch_splines.x)
set_tests_properties(test_batch_splines PROPERTIES
    LABELS "unit;performance"
    TIMEOUT 60
)

add_subdirectory(field_can)
add_subdirectory(magfie)
add_subdirectory(orbit)

##################################################
# ORBIT NETCDF OUTPUT SYSTEM TEST
##################################################

# Build test executable
add_executable(test_orbit_netcdf_output.x test_orbit_netcdf_output.f90)
target_link_libraries(test_orbit_netcdf_output.x simple)

# Create symlink to wout.nc in example directory
file(CREATE_LINK "${WOUT_FILE}" "${CMAKE_SOURCE_DIR}/examples/orbit_output_test/wout.nc" SYMBOLIC)

# System test: Run SIMPLE to generate orbits.nc, then verify
add_test(
    NAME orbit_netcdf_output
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/examples/orbit_output_test
        ${CMAKE_COMMAND} -E env
        $<TARGET_FILE:simple.x>
)
set_tests_properties(orbit_netcdf_output PROPERTIES
    TIMEOUT 120
    LABELS "system;integration"
)

# Verification test: Check that orbits.nc is valid
add_test(
    NAME orbit_netcdf_verify
    COMMAND test_orbit_netcdf_output.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(orbit_netcdf_verify PROPERTIES
    DEPENDS orbit_netcdf_output
    TIMEOUT 30
    LABELS "system;integration"
)

if(SIMPLE_ENABLE_PYTHON_TOOLS)
    add_test(
        NAME orbit_netcdf_plot
        COMMAND ${CMAKE_COMMAND} -E env python3 ${CMAKE_CURRENT_SOURCE_DIR}/plot_orbit.py 3
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(orbit_netcdf_plot PROPERTIES
        DEPENDS orbit_netcdf_output
        TIMEOUT 30
        LABELS "system;integration;plot;python"
    )
endif()

##################################################
# GOLDEN RECORD TESTS
##################################################

# Copy test cases to build directory for golden record tests
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record/boozer
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/golden_record/test_cases/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record/canonical
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/golden_record/test_cases/)

# Find all golden record test cases dynamically by scanning for simple.in files
file(GLOB GOLDEN_RECORD_CASES
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record
    ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record/*/simple.in
)

# Extract test case names (directory names)
set(GOLDEN_RECORD_TEST_CASES)
foreach(CASE_PATH ${GOLDEN_RECORD_CASES})
    get_filename_component(CASE_NAME ${CASE_PATH} DIRECTORY)
    list(APPEND GOLDEN_RECORD_TEST_CASES ${CASE_NAME})
endforeach()

if(GOLDEN_RECORD_TEST_CASES)
    message(STATUS "Golden record test cases found: ${GOLDEN_RECORD_TEST_CASES}")
else()
    message(WARNING "No golden record test cases found in test/golden_record/")
endif()

# Create golden record regression tests comparing against main branch
# These tests enforce strict numerical reproducibility - any change in output
# must be manually reviewed before merging to main
foreach(CASE ${GOLDEN_RECORD_TEST_CASES})
    add_test(
        NAME golden_record_${CASE}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record/golden_record.sh
            main ${CASE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(golden_record_${CASE} PROPERTIES
        TIMEOUT 300  # 5 minutes per case
        LABELS "golden_record;regression"
        ENVIRONMENT "GOLDEN_RECORD_BASE_DIR=${CMAKE_CURRENT_BINARY_DIR}/golden_record;GOLDEN_RECORD_RTOL=1.5e-7"
    )
endforeach()

# Golden record sanity test: tests the comparison system itself
# This ensures the golden record comparison scripts are working correctly
add_test(
    NAME golden_record_sanity
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record_sanity/test_golden_record_sanity_simple.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
)
set_tests_properties(golden_record_sanity PROPERTIES
    TIMEOUT 60  # 1 minute for sanity checks
    LABELS "golden_record;regression"
)
