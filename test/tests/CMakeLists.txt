
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/wout_unittests.nc)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/wout_unittests.nc
        ${CMAKE_CURRENT_BINARY_DIR}/wout.nc
        COPYONLY
    )
else()
    message(FATAL_ERROR "No .nc file, tests cannot be executed!")
endif()

add_pfunit_ctest(test_example
    TEST_SOURCES test_example.f90
    LINK_LIBRARIES SIMPLE_LIB)

add_pfunit_ctest(test_midlevel
    TEST_SOURCES test_midlevel.f90
    LINK_LIBRARIES SIMPLE_LIB)

add_pfunit_ctest(test_highlevel
    TEST_SOURCES test_highlevel.f90
    LINK_LIBRARIES SIMPLE_LIB)

add_pfunit_ctest(test_lowlevel
    TEST_SOURCES test_lowlevel.f90
    LINK_LIBRARIES SIMPLE_LIB)
    
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



