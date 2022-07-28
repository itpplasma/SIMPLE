
#ToDo rewrite the tests to match pFUnit
#ToDo this differently.
add_executable(test_boozer.x test_boozer.f90)
add_test(NAME test_boozer COMMAND test_boozer.x)

add_executable(test_coord_trans.x test_coord_trans.f90)
add_test(NAME test_coord_trans COMMAND test_coord_trans.x)

add_executable(test_magfie.x test_magfie.f90)
add_test(NAME test_magfie COMMAND test_magfie.x)

add_executable(test_new_can.x test_new_can.f90)
add_test(NAME test_magfie COMMAND test_magfie.x)

add_executable(test_orbits_vmec.x test_orbits_vmec.f90)
add_test(NAME test_orbits_vmec COMMAND test_orbits_vmec.x)

add_executable(test_parse_ants.x test_parse_ants.f90)
add_test(NAME test_parse_ants COMMAND test_parse_ants.x)

add_executable(test_poincare1.x test_poincare1.f90)
add_test(NAME test_poincare1 COMMAND test_poincare1.x)

add_executable(test_poiplot_classification.x test_poiplot_classification.f90)
add_test(NAME test_poiplot_classification COMMAND test_poiplot_classification.x)

add_executable(test_sympl.x test_sympl.f90)
add_test(NAME test_sympl COMMAND test_sympl.x)

add_executable(test_sympl_stell.x test_sympl_stell.f90)
add_test(NAME test_sympl_stell COMMAND test_sympl_stell.x)

add_executable(test_vmec.x test_vmec.f90)
add_test(NAME test_vmec COMMAND test_vmec.x)


