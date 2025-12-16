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
