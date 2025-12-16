add_executable (test_magfie_coils.x test_magfie_coils.f90)
target_link_libraries(test_magfie_coils.x simple)

add_executable (test_coils_field_compare.x test_coils_field_compare.f90)
target_link_libraries(test_coils_field_compare.x simple)

add_executable(test_magfie_refcoords_vmec.x test_magfie_refcoords_vmec.f90)
target_link_libraries(test_magfie_refcoords_vmec.x simple)

add_test(NAME test_magfie_refcoords_vmec
    COMMAND $<TARGET_FILE:test_magfie_refcoords_vmec.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
)
