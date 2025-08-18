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

# Python batch API tests for Issue #90
add_test(NAME test_batch_api
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_batch_api.py -v
)
set_tests_properties(test_batch_api PROPERTIES 
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}"
    LABELS "python;batch_api"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Performance and scalability tests (marked as slow)
add_test(NAME test_batch_api_performance
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_batch_api.py::TestPerformanceBenchmarking -v
)
set_tests_properties(test_batch_api_performance PROPERTIES 
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}"
    LABELS "python;batch_api;slow;performance"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Scalability tests (marked as slow and requiring large memory)
add_test(NAME test_batch_api_scalability
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_batch_api.py -k "scalability" -v
)
set_tests_properties(test_batch_api_scalability PROPERTIES 
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}"
    LABELS "python;batch_api;slow;scalability"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
