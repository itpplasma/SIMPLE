##################################################
# INTERNAL UNIT TEST STUFF
##################################################
add_library(testing_utils ${CMAKE_CURRENT_SOURCE_DIR}/utility.f90 ${CMAKE_CURRENT_SOURCE_DIR}/results.f90)

##################################################
# TEST EXECUTABLES
##################################################
add_executable(utility_test.x utility_test.f90)
target_link_libraries(utility_test.x testing_utils)

if (ENABLE_OPENMP)
    add_executable (test_sympl_tok.x test_sympl.f90)
    target_link_libraries(test_sympl_tok.x simple)
    add_executable (test_sympl.x test_sympl_stell.f90)
    target_link_libraries(test_sympl.x simple)
endif()

add_executable (test_coord_trans.x test_coord_trans.f90)
target_link_libraries(test_coord_trans.x simple)

add_executable (test_magfie.x test_magfie.f90)
target_link_libraries(test_magfie.x simple)
add_executable (test_poincare1.x test_poincare1.f90)
target_link_libraries(test_poincare1.x simple)
add_executable (test_poiplot_classification.x test_poiplot_classification.f90
)
target_link_libraries(test_poiplot_classification.x simple)
add_executable (test_parse_ants.x test_parse_ants.f90)
target_link_libraries(test_parse_ants.x simple)
add_executable (test_boozer.x test_boozer.f90)
target_link_libraries(test_boozer.x simple)
add_executable (test_orbits.x test_orbits_vmec.f90)
target_link_libraries(test_orbits.x simple)
add_executable (test_collis.x ../../src/test_collis.f90)
target_link_libraries(test_collis.x simple)

add_executable (test_samplers.x test_samplers.f90)
target_link_libraries(test_samplers.x simple)

add_executable (test_coordinates.x test_coordinates.f90)
target_link_libraries(test_coordinates.x simple)
add_test(NAME test_coordinates COMMAND test_coordinates.x vmec cyl)

add_subdirectory(field_can)
add_subdirectory(magfie)
add_subdirectory(orbit)
