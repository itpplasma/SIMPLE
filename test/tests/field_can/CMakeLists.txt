add_executable (test_field_can_meiss.x test_field_can_meiss.f90)
target_link_libraries(test_field_can_meiss.x simple)
add_test(NAME test_field_can_meiss_config
    COMMAND $<TARGET_FILE:test_field_can_meiss.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
set_tests_properties(test_field_can_meiss_config PROPERTIES
    LABELS "integration"
    RUN_SERIAL TRUE
    TIMEOUT 120)

add_executable (test_field_can_albert.x test_field_can_albert.f90)
target_link_libraries(test_field_can_albert.x simple)
add_test(NAME test_field_can_albert_config
    COMMAND $<TARGET_FILE:test_field_can_albert.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
set_tests_properties(test_field_can_albert_config PROPERTIES
    LABELS "integration"
    RUN_SERIAL TRUE
    TIMEOUT 120)

add_executable (test_field_can_albert_coils_diagnostic.x test_field_can_albert_coils_diagnostic.f90)
target_link_libraries(test_field_can_albert_coils_diagnostic.x simple)

add_executable (test_albert_intermediate.x test_albert_intermediate.f90)
target_link_libraries(test_albert_intermediate.x simple)

add_executable(test_coord_transform_roundtrip.x test_coord_transform_roundtrip.f90)
target_link_libraries(test_coord_transform_roundtrip.x simple)

add_test(NAME test_coord_transform_roundtrip
    COMMAND $<TARGET_FILE:test_coord_transform_roundtrip.x>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
set_tests_properties(test_coord_transform_roundtrip PROPERTIES
    LABELS "integration;slow"
    TIMEOUT 300)

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

# Per-volume Meiss canonical coordinates for SPECTRE (#439): construction identity
# on the committed libneo stellarator fixture, p_phi conservation on the committed
# axisymmetric two-volume tokamak fixture.
add_executable(test_spectre_sympl_volume.x test_spectre_sympl_volume.f90)
target_link_libraries(test_spectre_sympl_volume.x simple)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_construction)
add_test(NAME test_spectre_sympl_volume_construction
    COMMAND $<TARGET_FILE:test_spectre_sympl_volume.x> construction
        ${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_construction)
set_tests_properties(test_spectre_sympl_volume_construction PROPERTIES
    LABELS "integration"
    TIMEOUT 300)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_pphi)
add_test(NAME test_spectre_sympl_volume_pphi
    COMMAND $<TARGET_FILE:test_spectre_sympl_volume.x> pphi
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_pphi)
set_tests_properties(test_spectre_sympl_volume_pphi PROPERTIES
    LABELS "integration"
    TIMEOUT 300)

add_executable(test_spectre_multivolume_init.x test_spectre_multivolume_init.f90)
target_link_libraries(test_spectre_multivolume_init.x simple)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_multivolume_init)
add_test(NAME test_spectre_multivolume_init
    COMMAND $<TARGET_FILE:test_spectre_multivolume_init.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_multivolume_init)
set_tests_properties(test_spectre_multivolume_init PROPERTIES
    LABELS "integration"
    ENVIRONMENT "OMP_NUM_THREADS=1"
    TIMEOUT 60)

# Full orbit near SPECTRE interfaces (crossing_level = 3): GC <-> FO switching
# on the axisymmetric current-sheet fixture, with the momentum shell and the
# canonical toroidal momentum as witnesses.
add_executable(test_spectre_fo_interfaces.x test_spectre_fo_interfaces.f90)
target_link_libraries(test_spectre_fo_interfaces.x simple)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_fo_interfaces_rk45)
add_test(NAME test_spectre_fo_interfaces_rk45
    COMMAND $<TARGET_FILE:test_spectre_fo_interfaces.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_fo_interfaces_rk45)
set_tests_properties(test_spectre_fo_interfaces_rk45 PROPERTIES
    LABELS "integration;slow"
    TIMEOUT 600)

# Adversarial sharp-sheet events and gauge-consistent entry momentum.
add_executable(test_spectre_fo_events.x test_spectre_fo_events.f90)
target_link_libraries(test_spectre_fo_events.x simple)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_fo_events)
add_test(NAME test_spectre_fo_events
    COMMAND $<TARGET_FILE:test_spectre_fo_events.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_fo_events)
set_tests_properties(test_spectre_fo_events PROPERTIES
    LABELS "integration"
    TIMEOUT 60)

add_executable(test_spectre_fo_throughput.x test_spectre_fo_throughput.f90)
target_link_libraries(test_spectre_fo_throughput.x simple)
foreach(threads IN ITEMS 1 4)
    add_test(NAME test_spectre_fo_throughput_${threads}
        COMMAND $<TARGET_FILE:test_spectre_fo_throughput.x>
            ${CMAKE_CURRENT_SOURCE_DIR}/../data/ncsx_freeboundary.h5)
    set_tests_properties(test_spectre_fo_throughput_${threads} PROPERTIES
        LABELS "performance"
        ENVIRONMENT "OMP_NUM_THREADS=${threads}"
        RUN_SERIAL TRUE
        TIMEOUT 120)
endforeach()

add_executable(study_spectre_orbits.x
    ${CMAKE_SOURCE_DIR}/concept/spectre_loss_study/study_spectre_orbits.f90)
target_link_libraries(study_spectre_orbits.x simple)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_study)
add_test(NAME study_spectre_orbits
    COMMAND $<TARGET_FILE:study_spectre_orbits.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/ncsx_freeboundary.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_study)
set_tests_properties(study_spectre_orbits PROPERTIES
    LABELS "performance;research"
    RUN_SERIAL TRUE TIMEOUT 7200)

add_executable(test_spectre_fo_spline.x test_spectre_fo_spline.f90)
target_link_libraries(test_spectre_fo_spline.x simple)
foreach(fixture IN ITEMS tok2vol ncsx_freeboundary)
    add_test(NAME test_spectre_fo_spline_${fixture}
        COMMAND $<TARGET_FILE:test_spectre_fo_spline.x>
            ${CMAKE_CURRENT_SOURCE_DIR}/../data/${fixture}.h5)
    set_tests_properties(test_spectre_fo_spline_${fixture} PROPERTIES
        LABELS "integration"
        ENVIRONMENT "OMP_NUM_THREADS=4"
        RUN_SERIAL TRUE
        TIMEOUT 120)
endforeach()

add_test(NAME test_spectre_fo_spline_g3v3
    COMMAND $<TARGET_FILE:test_spectre_fo_spline.x>
        ${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5)
set_tests_properties(test_spectre_fo_spline_g3v3 PROPERTIES
    LABELS "integration" ENVIRONMENT "OMP_NUM_THREADS=4"
    RUN_SERIAL TRUE TIMEOUT 120)

# Diagnostic: zeroth-order GC -> FO -> GC switch round trip at NCSX interface
# points, without an FO step and over the exit gyrophase of a short full orbit;
# records d(mu)/mu, d(v_par), d(x_gc) and their gyroradius scaling. Asserts
# only completion and finite defects.
add_executable(test_spectre_switch_roundtrip.x test_spectre_switch_roundtrip.f90)
target_link_libraries(test_spectre_switch_roundtrip.x simple)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_switch_roundtrip)
add_test(NAME test_spectre_switch_roundtrip
    COMMAND $<TARGET_FILE:test_spectre_switch_roundtrip.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/ncsx_freeboundary.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_switch_roundtrip)
set_tests_properties(test_spectre_switch_roundtrip PROPERTIES
    LABELS "integration"
    TIMEOUT 300)

# Analytic magfie_spectre drift derivatives (bder, hcurl) must equal the central-
# difference algorithm they replaced to the FD truncation floor, at interior
# points of every volume of the committed tokamak and stellarator fixtures.
add_executable(test_magfie_spectre_analytic.x test_magfie_spectre_analytic.f90)
target_link_libraries(test_magfie_spectre_analytic.x simple)

add_test(NAME test_magfie_spectre_analytic_tok
    COMMAND $<TARGET_FILE:test_magfie_spectre_analytic.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_magfie_spectre_analytic_tok PROPERTIES
    LABELS "integration"
    TIMEOUT 120)

add_test(NAME test_magfie_spectre_analytic_stell
    COMMAND $<TARGET_FILE:test_magfie_spectre_analytic.x>
        ${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_magfie_spectre_analytic_stell PROPERTIES
    LABELS "integration"
    TIMEOUT 120)

add_executable(test_spectre_fo_gauss2.x test_spectre_fo_gauss2.f90)
target_link_libraries(test_spectre_fo_gauss2.x PRIVATE fo_gauss2_concept)
add_test(NAME test_spectre_fo_gauss2 COMMAND test_spectre_fo_gauss2.x
    ${CMAKE_SOURCE_DIR}/test/tests/data/tok2vol.h5 field)
set_tests_properties(test_spectre_fo_gauss2 PROPERTIES LABELS "integration" TIMEOUT 120)
# A non-exact interface gauge jump must be rejected during configuration.
add_test(NAME test_spectre_fo_gauss2_incompatible_gauge COMMAND test_spectre_fo_gauss2.x
    ${CMAKE_SOURCE_DIR}/test/tests/data/tok2vol.h5 incompatible)
set_tests_properties(test_spectre_fo_gauss2_incompatible_gauge PROPERTIES
    LABELS "unit" TIMEOUT 60
    PASS_REGULAR_EXPRESSION "interface normal fields disagree")
