

add_pfunit_ctest(test_boozer 
    TEST_SOURCES test_boozer.f90
	LINK_LIBRARIES SIMPLE
    LINK_LIBRARIES SIMPLE)

add_pfunit_ctest(test_coord_trans 
	TEST_SOURCES test_coord_trans.f90
	LINK_LIBRARIES SIMPLE)

add_pfunit_ctest(test_magfie 
	TEST_SOURCES test_magfie.f90
	LINK_LIBRARIES SIMPLE)

add_pfunit_ctest(test_magfie_unit
	TEST_SOURCES test_magfie.f90
	LINK_LIBRARIES SIMPLE)

add_pfunit_ctest(test_orbits_vmec 
	TEST_SOURCES test_orbits_vmec.f90
	LINK_LIBRARIES SIMPLE)

add_pfunit_ctest(test_parse_ants 
	TEST_SOURCES test_parse_ants.f90
	LINK_LIBRARIES SIMPLE)

add_pfunit_ctest(test_poincare1 
	TEST_SOURCES test_poincare1.f90
	LINK_LIBRARIES SIMPLE)

add_pfunit_ctest(test_poiplot_classification 
	TEST_SOURCES test_poiplot_classification.f90
	LINK_LIBRARIES SIMPLE)

add_pfunit_ctest(test_sympl 
	TEST_SOURCES test_sympl.f90
	LINK_LIBRARIES SIMPLE)

add_pfunit_ctest(test_sympl_stell 
	TEST_SOURCES test_sympl_stell.f90
	LINK_LIBRARIES SIMPLE)

add_pfunit_ctest(test_vmec 
	TEST_SOURCES test_vmec.f90
	LINK_LIBRARIES SIMPLE)


