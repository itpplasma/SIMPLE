##################################################
# INTERNAL UNIT TEST STUFF 
##################################################
add_compile_options(-cpp)
add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-cpp>)
add_library(testing_utils ${CMAKE_CURRENT_SOURCE_DIR}/utility.f90 ${CMAKE_CURRENT_SOURCE_DIR}/results.f90)

##################################################
# TEST EXECUTABLES
##################################################
add_executable(utility_test.x utility_test.f90)
target_link_libraries(utility_test.x testing_utils)



add_executable (test_coord_trans.x test_coord_trans.f90)
target_link_libraries(test_coord_trans.x SIMPLE_LIB)
add_executable (test_sympl_tok.x test_sympl.f90)
target_link_libraries(test_sympl_tok.x  SIMPLE_LIB)
add_executable (test_sympl.x test_sympl_stell.f90)
target_link_libraries(test_sympl.x  SIMPLE_LIB)
add_executable (test_magfie.x test_magfie.f90)
target_link_libraries(test_magfie.x  SIMPLE_LIB)
add_executable (test_poincare1.x test_poincare1.f90)
target_link_libraries(test_poincare1.x  SIMPLE_LIB)
add_executable (test_poiplot_classification.x test_poiplot_classification.f90
)
target_link_libraries(test_poiplot_classification.x  SIMPLE_LIB)
add_executable (alpha_lifetime_vmec_c.x
	../../SRC/binsrc.f90
	../../SRC/alpha_lifetime_c.f90
)
target_link_libraries(alpha_lifetime_vmec_c.x  SIMPLE_LIB)
add_executable (alpha_lifetime_old.x
	../../SRC/binsrc.f90
	../../SRC/alpha_lifetime.f90
)
target_link_libraries(alpha_lifetime_old.x  SIMPLE_LIB)
add_executable (test_parse_ants.x test_parse_ants.f90)
target_link_libraries(test_parse_ants.x  SIMPLE_LIB)
add_executable (test_boozer.x test_boozer.f90)
target_link_libraries(test_boozer.x  SIMPLE_LIB)
add_executable (test_orbits.x test_orbits_vmec.f90)
target_link_libraries(test_orbits.x  SIMPLE_LIB)
add_executable (test_collis.x ../../SRC/test_collis.f90)
target_link_libraries(test_collis.x  SIMPLE_LIB)



