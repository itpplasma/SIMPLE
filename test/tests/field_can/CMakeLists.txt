add_executable (test_field_can_meiss.x test_field_can_meiss.f90)
target_link_libraries(test_field_can_meiss.x simple)

add_executable (test_field_can_albert.x test_field_can_albert.f90)
target_link_libraries(test_field_can_albert.x simple)

add_executable (test_field_can_transforms.x test_field_can_transforms.f90)
target_link_libraries(test_field_can_transforms.x simple)
