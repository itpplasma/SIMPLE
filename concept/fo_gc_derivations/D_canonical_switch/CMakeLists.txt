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
