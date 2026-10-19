# Regenerate FortSym leaf files into the build tree and require them to match
# the checked-in sources. Arguments: GENERATOR, OUTPUTS and REFERENCES (lists
# in the generator's argument order).
list(GET OUTPUTS 0 first_output)
get_filename_component(output_dir "${first_output}" DIRECTORY)
file(MAKE_DIRECTORY "${output_dir}")
execute_process(COMMAND "${GENERATOR}" ${OUTPUTS} RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "Generator failed: ${GENERATOR}")
endif()
foreach(output reference IN ZIP_LISTS OUTPUTS REFERENCES)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${output}"
        "${reference}" RESULT_VARIABLE differ)
    if(NOT differ EQUAL 0)
        message(FATAL_ERROR "${reference} is stale; rebuild the regenerate_ target")
    endif()
endforeach()
