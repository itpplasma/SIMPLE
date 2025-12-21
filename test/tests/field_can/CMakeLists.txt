add_executable (test_field_can_meiss.x test_field_can_meiss.f90)
target_link_libraries(test_field_can_meiss.x simple)

add_executable (test_field_can_albert.x test_field_can_albert.f90)
target_link_libraries(test_field_can_albert.x simple)

add_executable (test_field_can_transforms.x test_field_can_transforms.f90)
target_link_libraries(test_field_can_transforms.x simple)

add_executable (test_field_can_albert_diagnostic.x test_field_can_albert_diagnostic.f90)
target_link_libraries(test_field_can_albert_diagnostic.x simple)

add_executable (test_field_can_albert_coils_diagnostic.x test_field_can_albert_coils_diagnostic.f90)
target_link_libraries(test_field_can_albert_coils_diagnostic.x simple)

add_executable (test_albert_intermediate.x test_albert_intermediate.f90)
target_link_libraries(test_albert_intermediate.x simple)

add_executable (test_field_can_diagnostics.x test_field_can_diagnostics.f90)
target_link_libraries(test_field_can_diagnostics.x simple)

add_executable(test_coord_transform_roundtrip.x test_coord_transform_roundtrip.f90)
target_link_libraries(test_coord_transform_roundtrip.x simple)

add_test(NAME test_coord_transform_roundtrip
    COMMAND $<TARGET_FILE:test_coord_transform_roundtrip.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
set_tests_properties(test_coord_transform_roundtrip PROPERTIES
    LABELS "integration;slow"
    TIMEOUT 300)

add_executable(test_albert_transform_diagnostic.x test_albert_transform_diagnostic.f90)
target_link_libraries(test_albert_transform_diagnostic.x simple)

add_executable(test_scaling_abstraction.x test_scaling_abstraction.f90)
target_link_libraries(test_scaling_abstraction.x simple)

add_test(NAME test_scaling_abstraction
    COMMAND $<TARGET_FILE:test_scaling_abstraction.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
set_tests_properties(test_scaling_abstraction PROPERTIES
    LABELS "integration;slow"
    TIMEOUT 300)

add_executable(test_boozer_spline_derivatives.x test_boozer_spline_derivatives.f90)
target_link_libraries(test_boozer_spline_derivatives.x simple)

add_test(NAME test_boozer_spline_derivatives
    COMMAND $<TARGET_FILE:test_boozer_spline_derivatives.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
set_tests_properties(test_boozer_spline_derivatives PROPERTIES
    LABELS "unit"
    TIMEOUT 180)
