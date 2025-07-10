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
