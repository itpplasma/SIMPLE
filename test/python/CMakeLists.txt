find_package(Python COMPONENTS Interpreter)

file(DOWNLOAD
    https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc
    ${CMAKE_CURRENT_BINARY_DIR}/wout.nc
)

# Disabled test_vmec - uncomment to re-enable
# add_test(NAME test_vmec
#     COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_vmec.py
# )
# set_tests_properties(test_vmec PROPERTIES 
#     ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}"
#     LABELS "slow"
# )

function(add_simple_api_test test_name test_node)
    add_test(NAME test_simple_api_${test_name}
        COMMAND ${Python_EXECUTABLE} -m pytest
            ${CMAKE_CURRENT_SOURCE_DIR}/test_simple_api.py::${test_node} -v
    )
    set_tests_properties(test_simple_api_${test_name} PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
        LABELS "python;simple_api;slow"
        TIMEOUT 120
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endfunction()

function(add_chartmap_api_test test_name test_node)
    add_test(NAME test_simple_api_${test_name}
        COMMAND ${Python_EXECUTABLE} -m pytest
            ${CMAKE_CURRENT_SOURCE_DIR}/test_simple_api.py::${test_node} -v
    )
    set_tests_properties(test_simple_api_${test_name} PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python;SIMPLE_TEST_CHARTMAP=${CMAKE_BINARY_DIR}/test/tests/test_boozer_chartmap.nc"
        LABELS "python;simple_api;chartmap;slow"
        TIMEOUT 120
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS "generate_test_boozer_chartmap_data"
    )
endfunction()

add_simple_api_test(integration_constants TestConstants::test_integration_constants)
add_simple_api_test(default_values TestConstants::test_default_values)
add_simple_api_test(init_basic TestInitialization::test_init_basic)
add_simple_api_test(init_parameters TestInitialization::test_init_with_parameters)
add_simple_api_test(init_repeated TestInitialization::test_multiple_reinitializations)
add_chartmap_api_test(init_chartmap TestInitialization::test_chartmap_init_selects_boozer_field)
add_simple_api_test(integrator_alias TestIntegratorSelection::test_alias_resolution)
add_simple_api_test(sampler_shape TestSampling::test_surface_sampler_shape)
add_simple_api_test(sampler_sizes TestSampling::test_sample_different_sizes)
add_simple_api_test(sampler_volume TestSampling::test_sample_volume)
add_simple_api_test(sample_surface TestSampleSurface::test_sample_surface)
add_simple_api_test(sample_copy TestSampleSurface::test_particles_can_be_copied)
add_simple_api_test(trace_results TestTraceOrbits::test_trace_returns_dict_results)
add_simple_api_test(trace_vmec TestTraceOrbits::test_trace_with_vmec_field)
add_simple_api_test(trace_boozer TestTraceOrbits::test_trace_with_boozer_field)
add_simple_api_test(trace_filters TestTraceOrbits::test_confined_and_lost_filters)
add_simple_api_test(trace_integrators TestTraceOrbits::test_different_integrators)
add_simple_api_test(trace_single TestTraceOrbits::test_trace_single_orbit)
add_simple_api_test(trace_trajectory TestTraceOrbits::test_trace_with_trajectory)
add_simple_api_test(canonical_frequencies
    TestTraceOrbits::test_canonical_frequencies)
add_simple_api_test(canonical_frequency_validation
    TestTraceOrbits::test_canonical_frequency_validation)
add_simple_api_test(trace_to_tip
    TestTraceOrbits::test_trace_to_tip)
add_simple_api_test(trace_to_cut_validation
    TestTraceOrbits::test_trace_to_cut_validation)
add_simple_api_test(trace_skipped TestTraceOrbits::test_trace_parallel_preserves_skipped_passing_sentinel)
add_simple_api_test(trace_early_loss TestTraceOrbits::test_trace_orbit_uses_backend_loss_time_for_early_loss)

add_test(NAME test_batch_api
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_batch_api.py -v
)
set_tests_properties(test_batch_api PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;batch_api;slow"
    TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_classification_api
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_classification_api.py -v
)
set_tests_properties(test_classification_api PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;classification;slow"
    TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_orbit_macro
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_orbit_macro.py -v
)
set_tests_properties(test_orbit_macro PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;orbit;slow"
    TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_examples
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_examples.py -v
)
set_tests_properties(test_examples PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;examples;slow"
    TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Performance and scalability tests (marked as slow)
add_test(NAME test_batch_api_performance
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_batch_api.py::TestPerformanceFramework::test_column_major_access_speed -v
)
set_tests_properties(test_batch_api_performance PROPERTIES 
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;batch_api;performance"
    TIMEOUT 120
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Scalability tests (marked as slow and requiring large memory)
add_test(NAME test_batch_api_scalability
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_batch_api.py::TestPerformanceFramework::test_column_major_access_speed -v
)
set_tests_properties(test_batch_api_scalability PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;batch_api;scalability"
    TIMEOUT 120
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# NetCDF results output tests
add_test(NAME test_netcdf_results
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_netcdf_results.py -v
)
set_tests_properties(test_netcdf_results PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;netcdf;slow"
    TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Chartmap boundary parameterization tests
add_test(NAME test_chartmap_boundary_param
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_chartmap_boundary_param.py -v
)
set_tests_properties(test_chartmap_boundary_param PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;chartmap;slow"
    TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Wall heat load engineering tests (no Fortran backend needed)
add_test(NAME test_engineering
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_engineering.py -v
)
set_tests_properties(test_engineering PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/python"
    LABELS "python;engineering"
    TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_coil_clearance
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_coil_clearance.py -v
)
set_tests_properties(test_coil_clearance PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/python"
    LABELS "python;engineering"
    TIMEOUT 120
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
