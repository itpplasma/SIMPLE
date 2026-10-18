file(READ "${CMAKE_CURRENT_LIST_DIR}/dependency_versions.json" _dependency_versions)
string(JSON _dependency_ref GET "${_dependency_versions}" "${DEPENDENCY}")
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E echo "${_dependency_ref}"
    COMMAND_ERROR_IS_FATAL ANY
)
