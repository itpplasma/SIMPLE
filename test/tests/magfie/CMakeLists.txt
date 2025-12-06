add_executable (test_magfie_coils.x test_magfie_coils.f90)
target_link_libraries(test_magfie_coils.x simple)

add_executable (test_coils_field_compare.x test_coils_field_compare.f90)
target_link_libraries(test_coils_field_compare.x simple)
