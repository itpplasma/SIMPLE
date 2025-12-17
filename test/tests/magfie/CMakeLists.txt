add_executable (test_magfie_coils.x test_magfie_coils.f90)
target_link_libraries(test_magfie_coils.x simple)

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

add_executable(test_orbit_chartmap_comparison.x test_orbit_chartmap_comparison.f90)
target_link_libraries(test_orbit_chartmap_comparison.x simple)
add_dependencies(test_orbit_chartmap_comparison.x chartmap_test_data)

add_executable(test_chartmap_refcoords_config.x test_chartmap_refcoords_config.f90)
target_link_libraries(test_chartmap_refcoords_config.x simple)
add_dependencies(test_chartmap_refcoords_config.x chartmap_test_data)

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

add_test(NAME test_hcurl_from_dhcov_splines
    COMMAND $<TARGET_FILE:test_hcurl_from_dhcov_splines.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)

add_test(NAME test_magfie_refcoords_splines
    COMMAND $<TARGET_FILE:test_magfie_refcoords_splines.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)

add_test(NAME test_orbit_refcoords_rk45
    COMMAND $<TARGET_FILE:test_orbit_refcoords_rk45.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)

add_test(NAME test_orbit_refcoords_rk45_plot
    COMMAND ${CMAKE_COMMAND} -E env python3
        ${CMAKE_CURRENT_SOURCE_DIR}/plot_orbit_refcoords_comparison.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
set_tests_properties(test_orbit_refcoords_rk45_plot PROPERTIES
    DEPENDS test_orbit_refcoords_rk45
    TIMEOUT 30
    LABELS "integration;plot"
)

add_test(NAME test_orbit_chartmap_comparison
    COMMAND $<TARGET_FILE:test_orbit_chartmap_comparison.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
set_tests_properties(test_orbit_chartmap_comparison PROPERTIES
    TIMEOUT 120
    LABELS "integration"
)

add_test(NAME test_chartmap_refcoords_config
    COMMAND $<TARGET_FILE:test_chartmap_refcoords_config.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
set_tests_properties(test_chartmap_refcoords_config PROPERTIES
    TIMEOUT 120
    LABELS "integration"
)

add_test(NAME test_orbit_chartmap_comparison_plot
    COMMAND ${CMAKE_COMMAND} -E env python3
        ${CMAKE_CURRENT_SOURCE_DIR}/plot_orbit_chartmap_comparison.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
set_tests_properties(test_orbit_chartmap_comparison_plot PROPERTIES
    DEPENDS test_orbit_chartmap_comparison
    TIMEOUT 60
    LABELS "integration;plot"
)

add_test(NAME test_orbit_chartmap_comparison_map2disc
    COMMAND $<TARGET_FILE:test_orbit_chartmap_comparison.x>
        wout.chartmap.map2disc.nc orbit_chartmap_comparison_map2disc.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
set_tests_properties(test_orbit_chartmap_comparison_map2disc PROPERTIES
    DEPENDS "generate_chartmap_map2disc"
    TIMEOUT 300
    LABELS "integration"
)

add_test(NAME test_orbit_chartmap_comparison_plot_map2disc
    COMMAND ${CMAKE_COMMAND} -E env python3
        ${CMAKE_CURRENT_SOURCE_DIR}/plot_orbit_chartmap_comparison.py
        orbit_chartmap_comparison_map2disc.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
set_tests_properties(test_orbit_chartmap_comparison_plot_map2disc PROPERTIES
    DEPENDS test_orbit_chartmap_comparison_map2disc
    TIMEOUT 60
    LABELS "integration;plot"
)

add_test(NAME test_chartmap_refcoords_config_map2disc
    COMMAND $<TARGET_FILE:test_chartmap_refcoords_config.x>
        wout.chartmap.map2disc.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
set_tests_properties(test_chartmap_refcoords_config_map2disc PROPERTIES
    DEPENDS "generate_chartmap_map2disc"
    TIMEOUT 120
    LABELS "integration"
)
