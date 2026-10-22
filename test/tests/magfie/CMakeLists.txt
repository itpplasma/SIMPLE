add_executable (test_magfie_coils.x test_magfie_coils.f90)
target_link_libraries(test_magfie_coils.x simple)
add_dependencies(test_magfie_coils.x ncsx_test_data)
add_test(NAME test_magfie_coils
    COMMAND $<TARGET_FILE:test_magfie_coils.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
set_tests_properties(test_magfie_coils PROPERTIES
    TIMEOUT 600
    LABELS "integration;slow"
)

add_executable (test_coils_field_compare.x test_coils_field_compare.f90)
target_link_libraries(test_coils_field_compare.x simple)

add_executable(test_magfie_refcoords_vmec.x test_magfie_refcoords_vmec.f90)
target_link_libraries(test_magfie_refcoords_vmec.x simple)

add_executable(test_jacobian_conventions_vmec.x test_jacobian_conventions_vmec.f90)
target_link_libraries(test_jacobian_conventions_vmec.x simple)

add_executable(test_hctr_from_metric_vmec.x test_hctr_from_metric_vmec.f90)
target_link_libraries(test_hctr_from_metric_vmec.x simple)

add_executable(test_hcurl_from_dhcov_splines.x test_hcurl_from_dhcov_splines.f90)
target_link_libraries(test_hcurl_from_dhcov_splines.x simple)

add_executable(test_magfie_refcoords_splines.x test_magfie_refcoords_splines.f90)
target_link_libraries(test_magfie_refcoords_splines.x simple)

add_executable(test_orbit_refcoords_rk45.x test_orbit_refcoords_rk45.f90)
target_link_libraries(test_orbit_refcoords_rk45.x simple)

add_executable(test_rk_pitch_zero_start.x test_rk_pitch_zero_start.f90)
target_link_libraries(test_rk_pitch_zero_start.x simple)

add_executable(test_orbit_chartmap_comparison.x test_orbit_chartmap_comparison.f90)
target_link_libraries(test_orbit_chartmap_comparison.x simple)
add_dependencies(test_orbit_chartmap_comparison.x chartmap_test_data)

add_executable(test_chartmap_wall_losses.x test_chartmap_wall_losses.f90)
target_link_libraries(test_chartmap_wall_losses.x simple)
add_dependencies(test_chartmap_wall_losses.x chartmap_test_data)

add_test(NAME test_magfie_refcoords_vmec
    COMMAND $<TARGET_FILE:test_magfie_refcoords_vmec.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)

add_test(NAME test_jacobian_conventions_vmec
    COMMAND $<TARGET_FILE:test_jacobian_conventions_vmec.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)

add_test(NAME test_hctr_from_metric_vmec
    COMMAND $<TARGET_FILE:test_hctr_from_metric_vmec.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)

# Sub-second VMEC-coordinate tests that exercise libneo directly; the libneo
# reverse-dependency gate runs this "smoke" set (ctest -L smoke).
set_property(TEST
    test_magfie_refcoords_vmec
    test_jacobian_conventions_vmec
    test_hctr_from_metric_vmec
    APPEND PROPERTY LABELS smoke)

add_test(NAME test_hcurl_from_dhcov_splines
    COMMAND $<TARGET_FILE:test_hcurl_from_dhcov_splines.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)

add_test(NAME test_magfie_refcoords_splines
    COMMAND $<TARGET_FILE:test_magfie_refcoords_splines.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)

# Peak RSS rounded up in GiB; these use the reduced coordinate fixture.
set_tests_properties(test_hcurl_from_dhcov_splines
    PROPERTIES RESOURCE_GROUPS "mem_gb:1")
set_tests_properties(test_magfie_refcoords_splines
    PROPERTIES RESOURCE_GROUPS "mem_gb:2")

add_test(NAME test_orbit_refcoords_rk45
    COMMAND $<TARGET_FILE:test_orbit_refcoords_rk45.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)

add_test(NAME test_rk_pitch_zero_start
    COMMAND $<TARGET_FILE:test_rk_pitch_zero_start.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)

add_test(NAME test_orbit_chartmap_comparison
    COMMAND $<TARGET_FILE:test_orbit_chartmap_comparison.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
set_tests_properties(test_orbit_chartmap_comparison PROPERTIES
    TIMEOUT 120
    LABELS "integration"
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set_tests_properties(test_orbit_chartmap_comparison PROPERTIES WILL_FAIL TRUE)
endif()

add_test(NAME test_chartmap_wall_losses
    COMMAND $<TARGET_FILE:test_chartmap_wall_losses.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
set_tests_properties(test_chartmap_wall_losses PROPERTIES
    TIMEOUT 300
    LABELS "integration;slow"
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set_tests_properties(test_chartmap_wall_losses PROPERTIES WILL_FAIL TRUE)
endif()

if(SIMPLE_ENABLE_PYTHON_TOOLS AND SIMPLE_HAVE_MAP2DISC)
    add_test(NAME test_orbit_chartmap_comparison_map2disc
        COMMAND $<TARGET_FILE:test_orbit_chartmap_comparison.x>
            wout.chartmap.map2disc.nc orbit_chartmap_comparison_map2disc.nc
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
    )
    set_tests_properties(test_orbit_chartmap_comparison_map2disc PROPERTIES
        DEPENDS "generate_chartmap_map2disc"
        TIMEOUT 300
        LABELS "integration;python"
        FIXTURES_REQUIRED chartmap_map2disc
    )

    add_test(NAME test_chartmap_wall_losses_map2disc
        COMMAND $<TARGET_FILE:test_chartmap_wall_losses.x>
            wout.chartmap.map2disc.nc chartmap_wall_losses_map2disc.nc
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
    )
    set_tests_properties(test_chartmap_wall_losses_map2disc PROPERTIES
        DEPENDS "generate_chartmap_map2disc"
        TIMEOUT 300
        LABELS "integration;python;slow"
        FIXTURES_REQUIRED chartmap_map2disc
    )
endif()
