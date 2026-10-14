execute_process(COMMAND "${GENERATOR}" "${OUTPUT}" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "Switch generator failed: ${status}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${OUTPUT}" "${EXPECTED}"
    RESULT_VARIABLE differs)
if(NOT differs EQUAL 0)
    message(FATAL_ERROR "Regenerated switch differs from committed source: ${OUTPUT}")
endif()
