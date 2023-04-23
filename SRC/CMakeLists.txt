

add_subdirectory(contrib)

set(SOURCES
    canonical_coordinates_mod.f90
    nctools_module.f90
    #    odeint_allroutines.f   # to use this, copy local version into SRC
    odeint_rkf45.f90       # interface to contrib/rkf45.f90
    magfie.f90
    boozer_converter.f90
    binsrc.f90
    plag_coeff.f90
    chamb_m.f90
    sub_alpha_lifetime_can.f90
    vmecinm_m.f90
    spline_vmec_data.f90
    spl_three_to_five.f90
    new_vmec_allocation_stuff.f90
    get_canonical_coordinates.f90
    testing.f90
    field_can.f90
    orbit_symplectic.f90
    orbit_symplectic_quasi.f90
    util.f90
    simple.f90
    bench.f90
    parse_ants.f90
    zzg.f90
    rng.f
    collis_alphas.f90
    params.f90
    sorting.f90
    check_orbit_type.f90
    find_bminmax.f90
)


add_library (SIMPLE_SRC SHARED ${SOURCES})
target_link_libraries(SIMPLE_SRC CONTRIB)

# Grab Python
find_package(Python REQUIRED
    COMPONENTS Interpreter Development NumPy)

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import numpy; print(numpy.__version__)"
    OUTPUT_VARIABLE NUMPY_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if (${NUMPY_VERSION} VERSION_GREATER_EQUAL "1.21.1")
	execute_process(
        COMMAND "${Python_EXECUTABLE}"
        -c "import numpy.f2py; print(numpy.f2py.get_include())"
        OUTPUT_VARIABLE F2PY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
	)
else()
    execute_process(
        COMMAND "${Python_EXECUTABLE}"
        -c "import numpy; print(numpy.get_include())"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(F2PY_INCLUDE_DIR "${NUMPY_INCLUDE_DIR}/../../f2py/src/")
endif()

include_directories(
    BEFORE
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
    ${F2PY_INCLUDE_DIR}
)

message(STATUS ${Python_INCLUDE_DIRS})
message(STATUS ${F2PY_INCLUDE_DIR})
message(STATUS ${Python_NumPy_INCLUDE_DIRS})
message(STATUS ${Python_LIBRARIES})

message(STATUS "Binary dir: ${CMAKE_BINARY_DIR}")

set(files_to_wrap
	simple.f90
	orbit_symplectic.f90
	field_can.f90
	get_canonical_coordinates.f90
	params.f90
	canonical_coordinates_mod.f90
	spline_vmec_data.f90
	bench.f90
)

set(preprocessed_sources
	${CMAKE_BINARY_DIR}/simple.f90.i
	${CMAKE_BINARY_DIR}/orbit_symplectic.f90.i
	${CMAKE_BINARY_DIR}/field_can.f90.i
	${CMAKE_BINARY_DIR}/get_canonical_coordinates.f90.i
	${CMAKE_BINARY_DIR}/params.f90.i
	${CMAKE_BINARY_DIR}/canonical_coordinates_mod.f90.i
	${CMAKE_BINARY_DIR}/spline_vmec_data.f90.i
	${CMAKE_BINARY_DIR}/bench.f90.i
)

add_custom_target(preprocessed_sources
	DEPENDS ${files_to_wrap}
)

foreach (file_to_wrap ${files_to_wrap})
	file(RELATIVE_PATH
		basename ${CMAKE_CURRENT_LIST_DIR}/ ${CMAKE_CURRENT_LIST_DIR}/${file_to_wrap}
		)
	add_custom_command(TARGET preprocessed_sources
		COMMAND gfortran -E -cpp
			-o ${CMAKE_BINARY_DIR}/${basename}.i
			${CMAKE_CURRENT_LIST_DIR}/${file_to_wrap}
		DEPENDS ${files_to_wrap}
	)
endforeach()


# Vars
set(f2py_module_name "_pysimple")
set(f2py_module_c "${f2py_module_name}module.c")

# Yes. f90.f90 due to preprocessor output filenames
set(fortran_src_files
	f90wrap_simple.f90.f90
	f90wrap_orbit_symplectic.f90.f90
	f90wrap_field_can.f90.f90
	f90wrap_get_canonical_coordinates.f90.f90
	f90wrap_params.f90.f90
	f90wrap_canonical_coordinates_mod.f90.f90
	f90wrap_bench.f90.f90
	f90wrap_toplevel.f90
)
set(generated_module_file "${f2py_module_name}.${Python_SOABI}")
set(f2py_wrappers "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}")


add_custom_command(
  OUTPUT ${f2py_wrappers}
  COMMAND f2py-f90wrap
                   ${fortran_src_files}
                   -m ${f2py_module_name}
                   --lower # Important
  DEPENDS ${fortran_src_files} # Fortran source
  )

add_custom_command(
	OUTPUT ${fortran_src_files}
	COMMAND f90wrap -m pysimple ${preprocessed_sources}
	DEPENDS preprocessed_sources
)

# Generate sources
add_custom_target(
  pysimple
  DEPENDS "${f2py_module_name}"
)

# # Set up target
add_library(${f2py_module_name} SHARED #EXCLUDE_FROM_ALL
   "${f2py_wrappers}"
   "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
   "${fortran_src_files}" # Fortran source(s)
)


# # Depend on sources
# add_dependencies(${f2py_module_name} genpyf)
target_link_libraries(${f2py_module_name} PRIVATE ${Python_LIBRARIES})

set_target_properties(
    ${f2py_module_name}
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "${generated_module_file}"
		SUFFIX ".so"
        LINKER_LANGUAGE C
    )

if (DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    add_custom_command(TARGET ${f2py_module_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_BINARY_DIR}/pysimple.py"
        "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/pysimple.py"
        COMMENT "Copying 'pysimple.py' file to '${CMAKE_LIBRARY_OUTPUT_DIRECTORY}'")
endif()

