

set(SOURCES
    test_boozer.f90
    test_coord_trans.f90
    test_magfie.f90
    test_new_can.f90
    test_orbits_vmec.90
    test_parse_ants.f90
    test_poincare1.f90
    test_poiplot_classification.f90
    test_sympl.f90
    test_sympl_stell.f90
    test_vmec.f90
)

add_library(TEST_CASES SHARED ${SOURCES})
