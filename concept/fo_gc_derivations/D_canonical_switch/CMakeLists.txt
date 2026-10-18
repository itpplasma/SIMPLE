# Consumer-owned generated transform. This library is linked only by concept
# diagnostics, never by SIMPLE's production orbit pipeline.
add_library(gc_switch_concept STATIC
    generated/gc_switch_generated.f90
    src/gc_fo_switch.f90
    src/spectre_jets.f90)
target_link_libraries(gc_switch_concept PUBLIC simple)

add_executable(study_canonical_switch.x tests/switch_study.f90)
target_link_libraries(study_canonical_switch.x PRIVATE gc_switch_concept)
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fields.json"
    CONTENT "{\"g3v3\":\"${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5\"}\n")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/smoke.in"
    "ncsx '${CMAKE_SOURCE_DIR}/test/tests/data/ncsx_freeboundary.h5' 1 0.5 0.4 0.3 0.3 0.7 0.0025 1 1000 20\n")
add_test(NAME study_canonical_switch
    COMMAND ${CMAKE_COMMAND}
        -DSTUDY=$<TARGET_FILE:study_canonical_switch.x>
        -DINPUT=${CMAKE_CURRENT_BINARY_DIR}/smoke.in
        -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/smoke.jsonl
        -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/native_study.cmake)
set_tests_properties(study_canonical_switch PROPERTIES
    LABELS "integration;research" TIMEOUT 1800
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

add_executable(test_switch_field_jets.x tests/test_switch_field_jets.f90)
target_link_libraries(test_switch_field_jets.x PRIVATE gc_switch_concept)
add_test(NAME test_switch_field_jets_ncsx COMMAND test_switch_field_jets.x
    "${CMAKE_SOURCE_DIR}/test/tests/data/ncsx_freeboundary.h5")
add_test(NAME test_switch_field_jets_g3v3 COMMAND test_switch_field_jets.x
    "${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5")
set_tests_properties(test_switch_field_jets_ncsx test_switch_field_jets_g3v3 PROPERTIES
    LABELS "integration;research" TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

add_executable(test_switch_pullback.x tests/test_switch_pullback.f90)
target_link_libraries(test_switch_pullback.x PRIVATE gc_switch_concept)
add_test(NAME test_switch_pullback_analytic COMMAND test_switch_pullback.x)
add_test(NAME test_switch_pullback_ncsx COMMAND test_switch_pullback.x
    "${CMAKE_SOURCE_DIR}/test/tests/data/ncsx_freeboundary.h5")
add_test(NAME test_switch_pullback_g3v3 COMMAND test_switch_pullback.x
    "${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5")
set_tests_properties(test_switch_pullback_analytic test_switch_pullback_ncsx
    test_switch_pullback_g3v3 PROPERTIES LABELS "integration;research" TIMEOUT 180
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

# Explicit research tooling only: ordinary builds require no CAS. Reuse the
# owning FortSym CMake build instead of compiling a private copy of its engine.
set(FORTSYM_BUILD_DIR "" CACHE PATH "FortSym CMake build for switch identity checks")
if(FORTSYM_BUILD_DIR)
    add_library(switch_fortsym SHARED IMPORTED)
    set_target_properties(switch_fortsym PROPERTIES
        IMPORTED_LOCATION "${FORTSYM_BUILD_DIR}/lib/libfortsym.so"
        INTERFACE_INCLUDE_DIRECTORIES "${FORTSYM_BUILD_DIR}/include")
    add_executable(check_switch_identities.x gen/check_switch_identities.f90)
    target_link_libraries(check_switch_identities.x PRIVATE switch_fortsym)
    add_executable(gen_switch1.x gen/gen_switch1.f90)
    target_link_libraries(gen_switch1.x PRIVATE switch_fortsym)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/results")
    add_test(NAME test_switch_symbolic_identities COMMAND check_switch_identities.x)
    set_tests_properties(test_switch_symbolic_identities PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" LABELS "research;symbolic"
        TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    add_test(NAME test_switch_regeneration
        COMMAND ${CMAKE_COMMAND} -DGENERATOR=$<TARGET_FILE:gen_switch1.x>
            -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/gc_switch_regenerated.f90
            -DEXPECTED=${CMAKE_CURRENT_SOURCE_DIR}/generated/gc_switch_generated.f90
            -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_regeneration.cmake)
    set_tests_properties(test_switch_regeneration PROPERTIES
        LABELS "research;symbolic" TIMEOUT 900
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
