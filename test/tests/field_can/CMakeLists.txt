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

# Symplectic multi-volume crossing (#441): long-run H and p_phi invariants of
# the composed flow-jump-flow map on the committed tok2vol fixture.
add_executable(test_spectre_sympl_crossing.x test_spectre_sympl_crossing.f90)
target_link_libraries(test_spectre_sympl_crossing.x simple)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_sympl_crossing)
add_test(NAME test_spectre_sympl_crossing_invariants
    COMMAND $<TARGET_FILE:test_spectre_sympl_crossing.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_sympl_crossing)
set_tests_properties(test_spectre_sympl_crossing_invariants PROPERTIES
    LABELS "integration"
    TIMEOUT 300)

# Full orbit near SPECTRE interfaces (crossing_level = 3): GC <-> FO switching
# on the axisymmetric current-sheet fixture, with the momentum shell and the
# canonical toroidal momentum as witnesses.
add_executable(test_spectre_fo_interfaces.x test_spectre_fo_interfaces.f90)
target_link_libraries(test_spectre_fo_interfaces.x simple)

foreach(fo_integrator boris rk45 gauss2)
    file(MAKE_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/spectre_fo_interfaces_${fo_integrator})
    add_test(NAME test_spectre_fo_interfaces_${fo_integrator}
        COMMAND $<TARGET_FILE:test_spectre_fo_interfaces.x>
            ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol.h5 ${fo_integrator}
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}/spectre_fo_interfaces_${fo_integrator})
    set_tests_properties(test_spectre_fo_interfaces_${fo_integrator} PROPERTIES
        LABELS "integration;slow"
        TIMEOUT 600)
    if(fo_integrator STREQUAL "gauss2")
        # Fixed fictitious steps spend most work resolving the layer. Allow
        # the long invariant run to finish on shared or low-core-count hosts.
        set_tests_properties(test_spectre_fo_interfaces_${fo_integrator}
            PROPERTIES TIMEOUT 1800)
    endif()
endforeach()

# Kicked Level-1 reflection (CA15): forbidden crossings relocate along the
# interface to the conjugate equal-B point at frozen (vpar, p, lambda, mu).
add_executable(test_spectre_kicked_reflection.x test_spectre_kicked_reflection.f90)
target_link_libraries(test_spectre_kicked_reflection.x simple)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_kicked_reflection)
add_test(NAME test_spectre_kicked_reflection
    COMMAND $<TARGET_FILE:test_spectre_kicked_reflection.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol_sheet.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_kicked_reflection)
set_tests_properties(test_spectre_kicked_reflection PROPERTIES
    LABELS "integration"
    TIMEOUT 300)

# Physical mu-scatter reference (ROADMAP section 4): full-orbit interface
# traversals against the conserved-mu GC sheet convention, two gyrostep scales.
add_executable(test_spectre_fo_mu_scatter.x test_spectre_fo_mu_scatter.f90)
target_link_libraries(test_spectre_fo_mu_scatter.x simple)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_fo_mu_scatter)
add_test(NAME test_spectre_fo_mu_scatter
    COMMAND $<TARGET_FILE:test_spectre_fo_mu_scatter.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol_sheet.h5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spectre_fo_mu_scatter)
set_tests_properties(test_spectre_fo_mu_scatter PROPERTIES
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

# Independent field/derivative oracles, analytic tangent map and stationary
# full-orbit energy envelope across repeated layer crossings.
add_executable(test_spectre_fo_gauss2.x test_spectre_fo_gauss2.f90)
target_link_libraries(test_spectre_fo_gauss2.x simple)
add_test(NAME test_spectre_fo_gauss2
    COMMAND $<TARGET_FILE:test_spectre_fo_gauss2.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol.h5)
set_tests_properties(test_spectre_fo_gauss2 PROPERTIES
    LABELS "integration;slow"
    TIMEOUT 600)

add_test(NAME test_spectre_fo_gauss2_ncsx
    COMMAND $<TARGET_FILE:test_spectre_fo_gauss2.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/ncsx_freeboundary.h5 field)
set_tests_properties(test_spectre_fo_gauss2_ncsx PROPERTIES
    LABELS "integration;slow"
    TIMEOUT 120)

add_test(NAME test_spectre_fo_gauss2_incompatible_gauge
    COMMAND $<TARGET_FILE:test_spectre_fo_gauss2.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../data/tok2vol.h5 incompatible)
set_tests_properties(test_spectre_fo_gauss2_incompatible_gauge PROPERTIES
    LABELS "integration"
    WILL_FAIL TRUE
    TIMEOUT 120)
