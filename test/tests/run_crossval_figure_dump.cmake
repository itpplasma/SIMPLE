# Run diag_crossval_figures.x and assert the analytic dump files exist with data.
# DIAG_EXE and OUT_DIR come from the ctest add_test -D flags. The VMEC files are
# optional (only written when wout.nc is in the working directory).

execute_process(
    COMMAND "${DIAG_EXE}" "${OUT_DIR}"
    RESULT_VARIABLE rc
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "diag_crossval_figures.x failed (rc=${rc}):\n${out}\n${err}")
endif()

# Trajectory dumps: many rows. The dt-sweep is one row per step size (5 + header).
set(trajectories
    analytic_cp.dat analytic_cpp_sym.dat analytic_cpp_var.dat
    analytic_gc.dat banana_pauli.dat banana_gc_drift.dat)
set(tables analytic_cpp_sym_dtsweep.dat)

foreach(f ${trajectories})
    set(path "${OUT_DIR}/${f}")
    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "missing dump file: ${path}")
    endif()
    file(STRINGS "${path}" lines)
    list(LENGTH lines n)
    if(n LESS 10)
        message(FATAL_ERROR "dump file too short (${n} lines): ${path}")
    endif()
endforeach()

foreach(f ${tables})
    set(path "${OUT_DIR}/${f}")
    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "missing dump file: ${path}")
    endif()
    file(STRINGS "${path}" lines)
    list(LENGTH lines n)
    if(n LESS 4)
        message(FATAL_ERROR "table dump too short (${n} lines): ${path}")
    endif()
endforeach()

message(STATUS "crossval figure dump produced ${required}")
