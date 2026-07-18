# Enforce no procedure-pointer trampolines in test executables. A pointer to an
# internal procedure builds a stack trampoline that faults on a non-executable
# stack; keep procedure-pointer targets (e.g. test magfie backends) as module
# procedures. Inherited by the field_can/magfie/orbit subdirectories below.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Wtrampolines>
        $<$<COMPILE_LANGUAGE:Fortran>:-Werror=trampolines>
    )
endif()

##################################################
# DOWNLOAD TEST DATA
##################################################
# Download VMEC test file to unified test_data directory (shared with Python tests)
set(TEST_DATA_DIR "${CMAKE_SOURCE_DIR}/test/test_data")
set(WOUT_FILE "${TEST_DATA_DIR}/wout.nc")
set(WOUT_URL "https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc")
set(WOUT_QH_FILE "${TEST_DATA_DIR}/wout_qh.nc")
set(WOUT_QH_URL "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc")

# Download if not present
if(NOT EXISTS "${WOUT_FILE}")
    file(MAKE_DIRECTORY "${TEST_DATA_DIR}")
    message(STATUS "Downloading VMEC test file to ${WOUT_FILE}...")
    file(DOWNLOAD "${WOUT_URL}" "${WOUT_FILE}")
    message(STATUS "Downloaded VMEC test file")
endif()

if(NOT EXISTS "${WOUT_QH_FILE}")
    file(MAKE_DIRECTORY "${TEST_DATA_DIR}")
    message(STATUS "Downloading QH VMEC file to ${WOUT_QH_FILE}...")
    file(DOWNLOAD "${WOUT_QH_URL}" "${WOUT_QH_FILE}")
endif()

# Create symlink in test binary directory for backward compatibility
file(CREATE_LINK "${WOUT_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/wout.nc" SYMBOLIC)
file(CREATE_LINK "${WOUT_QH_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/wout_qh.nc" SYMBOLIC)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

if(SIMPLE_DETERMINISTIC_FP AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_test(
        NAME test_deterministic_dependency_flags
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../python/test_deterministic_dependency_flags.py
            ${CMAKE_BINARY_DIR}/compile_commands.json
            ${libneo_SOURCE_DIR}
    )
    set_tests_properties(test_deterministic_dependency_flags PROPERTIES
        LABELS "unit;regression")
endif()

set(BOOZER_CHARTMAP_PYTHON "${Python3_EXECUTABLE}")
if(EXISTS "${CMAKE_SOURCE_DIR}/.venv/bin/python")
    set(BOOZER_CHARTMAP_PYTHON "${CMAKE_SOURCE_DIR}/.venv/bin/python")
endif()

set(SIMPLE_HAVE_GVEC OFF)
execute_process(
    COMMAND ${BOOZER_CHARTMAP_PYTHON} -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('gvec') else 1)"
    RESULT_VARIABLE GVEC_IMPORT_RESULT
)
if(GVEC_IMPORT_RESULT EQUAL 0)
    set(SIMPLE_HAVE_GVEC ON)
else()
    message(STATUS "gvec Python package not found; skipping GVEC chartmap QA test")
endif()

set(SIMPLE_HAVE_BOOZ_XFORM OFF)
execute_process(
    COMMAND ${BOOZER_CHARTMAP_PYTHON} -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('booz_xform') else 1)"
    RESULT_VARIABLE BOOZ_XFORM_IMPORT_RESULT
)
if(BOOZ_XFORM_IMPORT_RESULT EQUAL 0)
    set(SIMPLE_HAVE_BOOZ_XFORM ON)
else()
    message(STATUS "booz_xform Python package not found; skipping booz_xform chartmap test")
endif()

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
set(SIMPLE_HAVE_MAP2DISC OFF)
if(SIMPLE_ENABLE_PYTHON_TOOLS)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('map2disc') else 1)"
        RESULT_VARIABLE MAP2DISC_IMPORT_RESULT
    )
    if(MAP2DISC_IMPORT_RESULT EQUAL 0)
        set(SIMPLE_HAVE_MAP2DISC ON)
    else()
        message(STATUS "map2disc Python package not found; skipping map2disc-based chartmap tests")
    endif()
endif()

if(SIMPLE_ENABLE_PYTHON_TOOLS AND SIMPLE_HAVE_MAP2DISC)
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

    # libneo's Python package imports compiled extensions (e.g. _efit_to_boozer).
    # When building libneo as a subproject, those extensions are written to the
    # libneo binary directory.
    if(DEFINED libneo_BINARY_DIR)
        set(LIBNEO_PYTHON_BINARY_DIR "${libneo_BINARY_DIR}")
    else()
        set(LIBNEO_PYTHON_BINARY_DIR "${CMAKE_BINARY_DIR}/libneo")
    endif()

    add_test(NAME generate_chartmap_map2disc
        COMMAND ${CMAKE_COMMAND} -E env
            PYTHONPATH=${LIBNEO_PYTHON_DIR}:${LIBNEO_PYTHON_BINARY_DIR}
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
            PYTHONPATH=${LIBNEO_PYTHON_DIR}:${LIBNEO_PYTHON_BINARY_DIR}
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
            PYTHONPATH=${LIBNEO_PYTHON_DIR}:${LIBNEO_PYTHON_BINARY_DIR}
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
            PYTHONPATH=${LIBNEO_PYTHON_DIR}:${LIBNEO_PYTHON_BINARY_DIR}
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
            PYTHONPATH=${LIBNEO_PYTHON_DIR}:${LIBNEO_PYTHON_BINARY_DIR}
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
    # Stellarator symplectic + quasi-symplectic integrator driver. The quasi
    # branch is the only code path that exercises the root solver behind the
    # quasi steps (fortnum multiroot_hybrids on this branch, MINPACK hybrd1 on
    # main); registering it gives that path real golden coverage.
    add_test(NAME test_sympl_stell
        COMMAND test_sympl.x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(test_sympl_stell PROPERTIES
        LABELS "integration"
        RUN_SERIAL TRUE
        TIMEOUT 120)
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

    add_executable(test_orbit_frequencies.x test_orbit_frequencies.f90)
    target_link_libraries(test_orbit_frequencies.x simple)
    add_test(NAME test_orbit_frequencies
        COMMAND test_orbit_frequencies.x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(test_orbit_frequencies PROPERTIES TIMEOUT 120 LABELS "physics")

    # Axis-crossing regression (issue #370): near-axis trapped orbits must
    # never be lost at the axis and must conserve energy across crossings.
    add_executable(test_axis_crossing.x test_axis_crossing.f90)
    target_link_libraries(test_axis_crossing.x simple)
    add_test(NAME test_axis_crossing
        COMMAND test_axis_crossing.x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(test_axis_crossing PROPERTIES
        LABELS "integration"
        TIMEOUT 120)

    # Full-orbit (FO) Boris pusher: energy bound, near-axis crossing, no spurious
    # loss, on the reactor-scale Boozer field (wout.nc symlinked above).
    add_executable(test_fo_boris.x test_fo_boris.f90)
    target_link_libraries(test_fo_boris.x simple)
    add_test(NAME test_fo_boris
        COMMAND test_fo_boris.x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(test_fo_boris PROPERTIES
        LABELS "integration"
        TIMEOUT 120)
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

add_executable (test_profiles.x test_profiles.f90)
target_link_libraries(test_profiles.x simple)
add_test(NAME test_profiles COMMAND test_profiles.x)
set_tests_properties(test_profiles PROPERTIES LABELS "unit")

add_executable (test_samplers.x test_samplers.f90)
target_link_libraries(test_samplers.x simple)

add_executable(test_deterministic_seed.x test_deterministic_seed.f90)
target_link_libraries(test_deterministic_seed.x simple)
add_test(NAME test_deterministic_seed COMMAND test_deterministic_seed.x)
set_tests_properties(test_deterministic_seed PROPERTIES LABELS "unit")

add_executable(test_interface_crossing_ownership.x test_interface_crossing_ownership.f90)
target_link_libraries(test_interface_crossing_ownership.x simple)
add_test(NAME test_interface_crossing_ownership COMMAND test_interface_crossing_ownership.x)
set_tests_properties(test_interface_crossing_ownership PROPERTIES LABELS "unit")

add_executable(test_spectre_sheet_gc.x test_spectre_sheet_gc.f90)
target_link_libraries(test_spectre_sheet_gc.x simple)
add_test(NAME test_spectre_sheet_gc
    COMMAND test_spectre_sheet_gc.x
        ${CMAKE_CURRENT_SOURCE_DIR}/data/tok2vol.h5)
set_tests_properties(test_spectre_sheet_gc PROPERTIES LABELS "integration")

add_executable(test_spectre_fo_hybrid.x test_spectre_fo_hybrid.f90)
target_link_libraries(test_spectre_fo_hybrid.x simple)
add_test(NAME test_spectre_fo_hybrid
    COMMAND test_spectre_fo_hybrid.x
        ${CMAKE_CURRENT_SOURCE_DIR}/data/tok2vol.h5)
set_tests_properties(test_spectre_fo_hybrid PROPERTIES LABELS "integration")

add_executable(test_spectre_freeboundary.x test_spectre_freeboundary.f90)
target_link_libraries(test_spectre_freeboundary.x simple)
add_test(NAME test_spectre_freeboundary
    COMMAND test_spectre_freeboundary.x
        ${CMAKE_CURRENT_SOURCE_DIR}/data/ncsx_freeboundary.h5)
set_tests_properties(test_spectre_freeboundary PROPERTIES
    LABELS "integration"
    TIMEOUT 120)

add_executable(test_bminmax_lifecycle.x test_bminmax_lifecycle.f90)
target_link_libraries(test_bminmax_lifecycle.x simple testing_utils)
add_test(NAME test_bminmax_lifecycle
    COMMAND test_bminmax_lifecycle.x
        ${CMAKE_CURRENT_SOURCE_DIR}/bminmax_cache_cases.tsv)
set_tests_properties(test_bminmax_lifecycle PROPERTIES
    LABELS "unit"
    TIMEOUT 15)

add_executable(test_orbit_invariants.x test_orbit_invariants.f90)
target_link_libraries(test_orbit_invariants.x simple)
add_test(NAME test_orbit_invariants COMMAND test_orbit_invariants.x)
set_tests_properties(test_orbit_invariants PROPERTIES LABELS "unit")

if(SIMPLE_ENABLE_PYTHON_TOOLS)
    add_test(NAME test_bminmax_python
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_bminmax_python.py
            ${CMAKE_CURRENT_SOURCE_DIR}/bminmax_cache_cases.tsv
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(test_bminmax_python PROPERTIES
        LABELS "unit;python"
        TIMEOUT 15)

    add_test(NAME test_bminmax_driver
        COMMAND ${CMAKE_COMMAND} -E env
            CTEST_BINARY_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}
            ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_bminmax_driver.py
            $<TARGET_FILE:simple.x>
            ${WOUT_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(test_bminmax_driver PROPERTIES
        LABELS "integration;system"
        TIMEOUT 120)
endif()

add_executable (test_coordinates.x test_coordinates.f90)
target_link_libraries(test_coordinates.x simple)
add_test(NAME test_coordinates COMMAND test_coordinates.x vmec cyl)

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

    add_executable(test_chartmap_metadata.x test_chartmap_metadata.f90)
    target_link_libraries(test_chartmap_metadata.x simple testing_utils)
    add_test(NAME test_chartmap_metadata COMMAND test_chartmap_metadata.x)
    set_tests_properties(test_chartmap_metadata PROPERTIES LABELS "unit")

    # Chartmap reader metadata and A_phi radial-abscissa contract.
    add_executable(test_chartmap_aphi_abscissa.x test_chartmap_aphi_abscissa.f90)
    target_link_libraries(test_chartmap_aphi_abscissa.x simple)
    add_test(NAME test_chartmap_aphi_abscissa COMMAND test_chartmap_aphi_abscissa.x)
    set_tests_properties(test_chartmap_aphi_abscissa PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        LABELS "unit")

    # Boozer chartmap test: generate test file with a compiled NetCDF tool, then run Fortran test
    add_executable(generate_test_boozer_chartmap.x generate_test_boozer_chartmap.c)
    target_link_libraries(generate_test_boozer_chartmap.x netcdf m)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/test_boozer_chartmap.nc
        COMMAND $<TARGET_FILE:generate_test_boozer_chartmap.x> ${CMAKE_CURRENT_BINARY_DIR}/test_boozer_chartmap.nc
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS generate_test_boozer_chartmap.x
        COMMENT "Generating test Boozer chartmap NetCDF"
    )
    add_custom_target(generate_test_boozer_chartmap
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/test_boozer_chartmap.nc)
    add_test(NAME generate_test_boozer_chartmap_data
        COMMAND $<TARGET_FILE:generate_test_boozer_chartmap.x> ${CMAKE_CURRENT_BINARY_DIR}/test_boozer_chartmap.nc)
    set_tests_properties(generate_test_boozer_chartmap_data PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        LABELS "unit")
    add_executable(test_boozer_chartmap.x test_boozer_chartmap.f90)
    target_link_libraries(test_boozer_chartmap.x simple)
    add_dependencies(test_boozer_chartmap.x generate_test_boozer_chartmap)
    add_test(NAME test_boozer_chartmap COMMAND test_boozer_chartmap.x)
    set_tests_properties(test_boozer_chartmap PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS "generate_test_boozer_chartmap_data"
        LABELS "unit")
    add_test(NAME test_chartmap_startmode1
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_chartmap_startmode1.py
            $<TARGET_FILE:simple.x>
            ${CMAKE_CURRENT_BINARY_DIR}/test_boozer_chartmap.nc)
    set_tests_properties(test_chartmap_startmode1 PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        LABELS "integration;python"
        DEPENDS "generate_test_boozer_chartmap_data")

    add_executable(test_chartmap_scaling.x test_chartmap_scaling.f90)
    target_link_libraries(test_chartmap_scaling.x simple)
    add_dependencies(test_chartmap_scaling.x generate_test_boozer_chartmap)
    add_test(NAME test_chartmap_scaling COMMAND test_chartmap_scaling.x)
    set_tests_properties(test_chartmap_scaling PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS "generate_test_boozer_chartmap_data"
        LABELS "unit")

    # Boozer chartmap roundtrip test: requires wout.nc
    add_executable(test_boozer_chartmap_roundtrip.x test_boozer_chartmap_roundtrip.f90)
    target_link_libraries(test_boozer_chartmap_roundtrip.x simple)
    add_test(NAME test_boozer_chartmap_roundtrip
        COMMAND test_boozer_chartmap_roundtrip.x
            ${WOUT_FILE}
            ${CMAKE_CURRENT_BINARY_DIR}/roundtrip_test.nc
            export
            ${CMAKE_CURRENT_BINARY_DIR}/boozer_chartmap_roundtrip
            8.0e-5
            2.0e-9)
    set_tests_properties(test_boozer_chartmap_roundtrip PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        LABELS "slow")

    if(SIMPLE_HAVE_GVEC)
        add_test(NAME test_boozer_chartmap_gvec_qa
            COMMAND ${BOOZER_CHARTMAP_PYTHON}
                ${CMAKE_CURRENT_SOURCE_DIR}/run_gvec_qa_roundtrip.py)
        set_tests_properties(test_boozer_chartmap_gvec_qa PROPERTIES
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            LABELS "slow;python"
            TIMEOUT 7200)
    endif()

    add_test(NAME test_gvec_chartmap_converter_static
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_gvec_chartmap_converter_static.py)
    set_tests_properties(test_gvec_chartmap_converter_static PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        LABELS "unit;python")

    # Export tool for e2e test
    add_executable(export_boozer_chartmap_tool.x export_boozer_chartmap_tool.f90)
    target_link_libraries(export_boozer_chartmap_tool.x simple)

    # Axis-loss regression (issue #370): zero symplectic-only losses on
    # identical near-axis starts, symplectic Euler vs RK Boozer reference.
    add_test(NAME test_axis_loss_consistency
        COMMAND ${BOOZER_CHARTMAP_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_axis_loss_consistency.py)
    set_tests_properties(test_axis_loss_consistency PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        LABELS "system;python"
        TIMEOUT 1800)

    # E2E test: VMEC-Boozer vs chartmap confined fractions
    add_test(NAME test_e2e_boozer_chartmap
        COMMAND ${BOOZER_CHARTMAP_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_e2e_boozer_chartmap.py)
    set_tests_properties(test_e2e_boozer_chartmap PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        LABELS "slow;plot;python"
        ENVIRONMENT "SIMPLE_BOOZER_CASES=QA"
        TIMEOUT 7200)

    # booz_xform -> Boozer chartmap converter: spectral reference, golden
    # cross-check against export_boozer_chartmap, end-to-end orbit comparison.
    if(SIMPLE_HAVE_BOOZ_XFORM)
        add_test(NAME test_booz_xform_chartmap
            COMMAND ${BOOZER_CHARTMAP_PYTHON}
                ${CMAKE_CURRENT_SOURCE_DIR}/test_booz_xform_chartmap.py)
        set_tests_properties(test_booz_xform_chartmap PROPERTIES
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            LABELS "slow;python"
            TIMEOUT 7200)
    endif()


    if(SIMPLE_ENABLE_CGAL)
        add_executable(test_stl_wall_intersection.x test_stl_wall_intersection.f90)
        target_link_libraries(test_stl_wall_intersection.x simple testing_utils)
        add_test(NAME test_stl_wall_intersection COMMAND test_stl_wall_intersection.x)
        set_tests_properties(test_stl_wall_intersection PROPERTIES LABELS "unit;wall")
    endif()

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

add_executable(test_newton_solver_status.x test_newton_solver_status.f90)
target_link_libraries(test_newton_solver_status.x simple)
add_test(NAME test_newton_solver_status COMMAND test_newton_solver_status.x)

add_executable(test_boundary_event_config.x test_boundary_event_config.f90)
target_link_libraries(test_boundary_event_config.x simple)
add_test(NAME test_boundary_event_config_valid
    COMMAND test_boundary_event_config.x valid)
add_test(NAME test_boundary_event_config_fraction_nan
    COMMAND test_boundary_event_config.x fraction_nan)
add_test(NAME test_boundary_event_config_fraction_inf
    COMMAND test_boundary_event_config.x fraction_inf)
add_test(NAME test_boundary_event_config_radial_nan
    COMMAND test_boundary_event_config.x radial_nan)
set_tests_properties(
    test_boundary_event_config_valid
    test_boundary_event_config_fraction_nan
    test_boundary_event_config_fraction_inf
    test_boundary_event_config_radial_nan
    PROPERTIES LABELS "unit")
set_tests_properties(
    test_boundary_event_config_fraction_nan
    test_boundary_event_config_fraction_inf
    test_boundary_event_config_radial_nan
    PROPERTIES WILL_FAIL TRUE)

add_executable(test_field_base.x test_field_base.f90)
target_link_libraries(test_field_base.x simple)
add_test(NAME test_field_base COMMAND test_field_base.x)

add_executable(test_ntcut_overflow.x test_ntcut_overflow.f90)
target_link_libraries(test_ntcut_overflow.x simple)
add_test(NAME test_ntcut_overflow COMMAND test_ntcut_overflow.x)

# Sub-second base-interface tests for the libneo reverse-dependency gate
# (ctest -L smoke).
set_property(TEST
    test_lapack_interfaces
    test_orbit_symplectic_base
    test_field_base
    test_ntcut_overflow
    APPEND PROPERTY LABELS smoke)

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
    RUN_SERIAL TRUE
    TIMEOUT 120)

# Chartmap R,Z consistency tests
add_executable(test_chartmap_rz_consistency.x test_chartmap_rz_consistency.f90)
target_link_libraries(test_chartmap_rz_consistency.x simple)
add_dependencies(test_chartmap_rz_consistency.x ncsx_test_data)

add_test(NAME test_chartmap_rz_consistency
    COMMAND test_chartmap_rz_consistency.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_chartmap_rz_consistency PROPERTIES
    LABELS "unit"
    RUN_SERIAL TRUE
    TIMEOUT 120)

# Batch spline optimization test
add_executable(test_batch_splines.x ../test_batch_splines.f90)
target_link_libraries(test_batch_splines.x simple testing_utils)
add_test(NAME test_batch_splines COMMAND test_batch_splines.x)
set_tests_properties(test_batch_splines PROPERTIES
    LABELS "unit;performance"
    TIMEOUT 60
)

# GPU device-evaluation test: batch spline eval inside an OpenACC kernel must
# match the host result. Only meaningful when OpenACC offload is enabled.
if(SIMPLE_ENABLE_OPENACC)
    add_executable(test_batch_splines_device.x test_batch_splines_device.f90)
    target_link_libraries(test_batch_splines_device.x simple testing_utils)
    add_test(NAME test_batch_splines_device COMMAND test_batch_splines_device.x)
    set_tests_properties(test_batch_splines_device PROPERTIES
        LABELS "unit;gpu"
        TIMEOUT 120
    )
endif()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gpu_orbit_bench)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_gpu_orbit_bench.in
    ${CMAKE_CURRENT_BINARY_DIR}/gpu_orbit_bench/test_gpu_orbit_bench.in
    COPYONLY)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_gpu_orbit_bench.start.txt
    ${CMAKE_CURRENT_BINARY_DIR}/gpu_orbit_bench/start.dat
    COPYONLY)

add_test(NAME test_gpu_orbit_bench
    COMMAND ${CMAKE_COMMAND} -E env
        SIMPLE_GPU_BENCH=1
        $<TARGET_FILE:simple.x>
        test_gpu_orbit_bench.in
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_gpu_orbit_bench PROPERTIES
    LABELS "unit"
    PASS_REGULAR_EXPRESSION "loss-step mismatches = *0"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gpu_orbit_bench
    TIMEOUT 120)

# Macrostep log grid tests (issue: performance regression with log grid)
add_executable(test_macrostep_log_grid.x test_macrostep_log_grid.f90)
target_link_libraries(test_macrostep_log_grid.x simple)
add_test(NAME test_macrostep_log_grid COMMAND test_macrostep_log_grid.x)
set_tests_properties(test_macrostep_log_grid PROPERTIES
    LABELS "unit;performance"
    TIMEOUT 60
)

# Restart I/O unit test
add_executable(test_restart_io.x test_restart_io.f90)
target_link_libraries(test_restart_io.x simple)
add_test(NAME test_restart_io
    COMMAND test_restart_io.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_restart_io PROPERTIES
    LABELS "unit"
    RESOURCE_LOCK simple_result_files
    TIMEOUT 15)

add_executable(test_result_accounting.x test_result_accounting.f90)
target_link_libraries(test_result_accounting.x simple)
add_test(NAME test_result_accounting
    COMMAND test_result_accounting.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_result_accounting PROPERTIES
    LABELS "unit"
    RESOURCE_LOCK simple_result_files
    TIMEOUT 15)

# Restart integration test (TEST field, deterministic)
add_test(NAME test_restart_integration
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_restart.py
        $<TARGET_FILE:simple.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_restart_integration PROPERTIES
    LABELS "integration;python"
    TIMEOUT 60)

# SPECTRE symplectic field-line Poincare app vs SPECTRE reference sections (#437).
# The fixture is the committed libneo G3V3L3Fi HDF5; do not copy it into SIMPLE.
add_test(NAME test_spectre_poincare
    COMMAND ${BOOZER_CHARTMAP_PYTHON}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_spectre_poincare.py
        $<TARGET_FILE:spectre_poincare.x>
        ${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5
        ${CMAKE_CURRENT_SOURCE_DIR}/data/poincare.npz
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_spectre_poincare PROPERTIES
    LABELS "integration;python"
    TIMEOUT 120)

add_executable(test_field_line_midpoint.x test_field_line_midpoint.f90)
target_link_libraries(test_field_line_midpoint.x simple)
add_test(NAME test_field_line_midpoint
    COMMAND test_field_line_midpoint.x
        ${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_field_line_midpoint PROPERTIES
    LABELS "unit;spectre"
    TIMEOUT 15)

# SPECTRE guiding-center RK45 tracing per volume with boundary-stop (#438).
# Uses simple.x with the committed libneo SPECTRE HDF5 fixture; VMEC-free.
add_test(NAME test_spectre_gc_rk45
    COMMAND ${BOOZER_CHARTMAP_PYTHON}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_spectre_gc_rk45.py
        $<TARGET_FILE:simple.x>
        ${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_spectre_gc_rk45 PROPERTIES
    LABELS "integration;python"
    TIMEOUT 120)

add_test(NAME test_spectre_startmode2
    COMMAND ${BOOZER_CHARTMAP_PYTHON}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_spectre_startmode2.py
        $<TARGET_FILE:simple.x>
        ${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_spectre_startmode2 PROPERTIES
    LABELS "integration;python"
    TIMEOUT 120)

# SPECTRE per-volume symplectic guiding-center tracing (#439): bounded non-secular
# energy, cross-scheme consistency, and RK45 cross-check on the committed fixture.
add_test(NAME test_spectre_sympl_volume
    COMMAND ${BOOZER_CHARTMAP_PYTHON}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_spectre_sympl_volume.py
        $<TARGET_FILE:simple.x>
        ${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_spectre_sympl_volume PROPERTIES
    LABELS "integration;python"
    RUN_SERIAL TRUE
    TIMEOUT 300)

# SPECTRE Level-0 interface crossing map with event location (#443). tok2vol is
# committed alongside the test; G3V3L3Fi is the committed libneo fixture.
add_test(NAME test_spectre_crossing_l0
    COMMAND ${BOOZER_CHARTMAP_PYTHON}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_spectre_crossing_l0.py
        $<TARGET_FILE:simple.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/data/tok2vol.h5
        ${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_spectre_crossing_l0 PROPERTIES
    LABELS "integration;python"
    TIMEOUT 300)

# SPECTRE Level-1 interface crossing refraction map (#440). Runs tok2vol under
# both crossing_level = 1 (default) and 0 and checks the map from the crossing
# log (generator components and lambda_k).
add_test(NAME test_spectre_crossing_l1
    COMMAND ${BOOZER_CHARTMAP_PYTHON}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_spectre_crossing_l1.py
        $<TARGET_FILE:simple.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/data/tok2vol.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_spectre_crossing_l1 PROPERTIES
    LABELS "integration;python"
    TIMEOUT 300)

# Symplectic multi-volume interface crossing (#441): exact-landing accuracy,
# reflections, RK45 cross-path loss consistency, and step-halving order on the
# committed tok2vol fixture. The long-run invariant witness runs separately as
# test_spectre_sympl_crossing_invariants (field_can/).
add_test(NAME test_spectre_sympl_crossing
    COMMAND ${BOOZER_CHARTMAP_PYTHON}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_spectre_sympl_crossing.py
        $<TARGET_FILE:simple.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/data/tok2vol.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_spectre_sympl_crossing PROPERTIES
    LABELS "integration;python"
    RUN_SERIAL TRUE
    TIMEOUT 300)

# SPECTRE construction-grid memory knob and axisymmetric phi auto-clamp (#442):
# the clamp triggers on axisymmetric tok2vol, peak RSS falls below the unclamped
# phi=32 build, and a coarse grid still reproduces the per-volume canonical chart.
add_test(NAME test_spectre_construction_grid
    COMMAND ${BOOZER_CHARTMAP_PYTHON}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_spectre_construction_grid.py
        $<TARGET_FILE:simple.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/data/tok2vol.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_spectre_construction_grid PROPERTIES
    LABELS "integration;python"
    TIMEOUT 120)

# SPECTRE end-to-end validation (#442): losses and accounting for both crossing
# maps, trapped/passing mirror-criterion split, guiding-center crossing mu
# conservation, and fixed-seed determinism on the committed tok2vol.
add_test(NAME test_spectre_validation
    COMMAND ${BOOZER_CHARTMAP_PYTHON}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_spectre_validation.py
        $<TARGET_FILE:simple.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/data/tok2vol.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_spectre_validation PROPERTIES
    LABELS "integration;python"
    RUN_SERIAL TRUE
    TIMEOUT 300)

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
# NCSX GC Cash-Karp COMPARISON TEST (VMEC vs MEISS)
##################################################

add_test(
    NAME gc_ck_ncsx_vmec_vs_meiss
    COMMAND ${CMAKE_COMMAND} -E env python3
        ${CMAKE_SOURCE_DIR}/examples/compare_ncsx_meiss_vmec_ck.py
        ${CMAKE_BINARY_DIR}/test/ncsx_meiss_vmec_ck
)
set_tests_properties(gc_ck_ncsx_vmec_vs_meiss PROPERTIES
    TIMEOUT 600
    LABELS "system;integration;plot;slow"
)

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
# must be manually reviewed before merging to main.
#
# The GOLDEN_RECORD_SKIP_CASES env var is still honored by
# compare_golden_results.sh for future intentional-divergence bug fixes;
# pass it through ENVIRONMENT when needed.
foreach(CASE ${GOLDEN_RECORD_TEST_CASES})
    add_test(
        NAME golden_record_${CASE}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../golden_record/golden_record.sh
            main ${CASE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(golden_record_${CASE} PROPERTIES
        TIMEOUT 1200  # Includes initial main clone/build on CI
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

# Chartmap Meiss debug test (rho_convention, scaling)
add_executable(test_chartmap_meiss_debug.x test_chartmap_meiss_debug.f90)
target_link_libraries(test_chartmap_meiss_debug.x simple)

add_test(NAME test_chartmap_meiss_debug
    COMMAND test_chartmap_meiss_debug.x ${CMAKE_CURRENT_BINARY_DIR}/wout.chartmap.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_chartmap_meiss_debug PROPERTIES
    LABELS "unit;chartmap"
    TIMEOUT 30
    FIXTURES_REQUIRED vmec_test_data)
