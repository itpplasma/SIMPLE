# Run simple.x on one and on two MPI ranks in separate directories and require
# that rank 0 traced only its share and wrote byte-identical result files.
# Expects MPIEXEC, NPFLAG, SIMPLE, INPUT, WOUT and WORK.
execute_process(COMMAND ${MPIEXEC} --version OUTPUT_VARIABLE version
    ERROR_VARIABLE version)
set(preflags "")
if(version MATCHES "Open MPI|OpenRTE")
    # A CI allocation may grant fewer slots than ranks.
    set(preflags --oversubscribe)
endif()

foreach(np 1 2)
    set(dir ${WORK}/np${np})
    file(REMOVE_RECURSE ${dir})
    file(MAKE_DIRECTORY ${dir})
    configure_file(${INPUT} ${dir}/simple.in COPYONLY)
    file(CREATE_LINK ${WOUT} ${dir}/wout.nc SYMBOLIC)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=2
            ${MPIEXEC} ${NPFLAG} ${np} ${preflags} ${SIMPLE} simple.in
        WORKING_DIRECTORY ${dir}
        RESULT_VARIABLE rc
        OUTPUT_VARIABLE out_${np}
        ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "simple.x on ${np} ranks failed (${rc}):\n${out_${np}}\n${err}")
    endif()
endforeach()

if(NOT out_2 MATCHES "MPI ranks = 2, particles on this rank = 8")
    message(FATAL_ERROR "two ranks did not split the 16 particles:\n${out_2}")
endif()

foreach(f start.dat times_lost.dat orbit_exit_code.dat confined_fraction.dat
          unresolved_fraction.dat avg_inverse_t_lost.dat)
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
        ${WORK}/np1/${f} ${WORK}/np2/${f} RESULT_VARIABLE differs)
    if(NOT differs EQUAL 0)
        message(FATAL_ERROR "${f} differs between one and two ranks")
    endif()
endforeach()
message(STATUS "one and two ranks wrote identical results")
