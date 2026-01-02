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

# Clean Python API tests
add_test(NAME test_simple_api
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_simple_api.py -v
)
set_tests_properties(test_simple_api PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;simple_api;slow"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_batch_api
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_batch_api.py -v
)
set_tests_properties(test_batch_api PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;batch_api;slow"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_classification_api
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_classification_api.py -v
)
set_tests_properties(test_classification_api PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;classification;slow"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_orbit_macro
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_orbit_macro.py -v
)
set_tests_properties(test_orbit_macro PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;orbit;slow"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_examples
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_examples.py -v
)
set_tests_properties(test_examples PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;examples;slow"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Performance and scalability tests (marked as slow)
add_test(NAME test_batch_api_performance
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_batch_api.py::TestPerformanceFramework::test_column_major_access_speed -v
)
set_tests_properties(test_batch_api_performance PROPERTIES 
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;batch_api;slow;performance"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Scalability tests (marked as slow and requiring large memory)
add_test(NAME test_batch_api_scalability
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_batch_api.py::TestPerformanceFramework::test_column_major_access_speed -v
)
set_tests_properties(test_batch_api_scalability PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;batch_api;slow;scalability"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# NetCDF results output tests
add_test(NAME test_netcdf_results
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_netcdf_results.py -v
)
set_tests_properties(test_netcdf_results PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python"
    LABELS "python;netcdf;slow"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
