# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.7)
project (SIMPLE Fortran C)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
#    set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_MACOSX_RPATH 1)

find_program(NF_CONFIG "nc-config")

if (NF_CONFIG)
execute_process(COMMAND "nc-config" --prefix
		OUTPUT_VARIABLE NFPREFIX)
else()
message(SEND_ERROR "nc-config not found. Please install libnetcdff-dev")
endif()

string(STRIP ${NFPREFIX} NFPREFIX)
set(NFINC ${NFPREFIX}/include)
set(NFLIBS ${NFPREFIX}/lib)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)


message(STATUS "CMake build type: " ${CMAKE_BUILD_TYPE})
message(STATUS "NetCDF include path: " ${NFINC})
message(STATUS "NetCDF lib path: " ${NFLIBS})
message(STATUS "MPI include path: " ${MPI_Fortran_INCLUDE_DIRS})

include_directories ($ENV{NETCDFF_INCLUDE} ${NFINC} ${MPI_Fortran_INCLUDE_DIRS})
link_directories ($ENV{NETCDF_LIB} $ENV{NETCDFF_LIB} ${NFLIBS} $ENV{HOME}/.local/lib)

add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-cpp>)
add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fopenmp>)
#add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008>)
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	add_compile_options(-mkl -qopenmp -warn all,nounused -stand f2008)
	link_libraries("-mkl -qopenmp")
	if (CMAKE_BUILD_TYPE MATCHES Debug)
		add_compile_options(-O0 -g -traceback)
	elseif (CMAKE_BUILD_TYPE MATCHES Profile)
		add_compile_options(-O2 -g -shared-intel -debug inline-debug-info
			-D TBB_USE_THREADING_TOOLS -qopenmp-link dynamic -parallel-source-info=2)
	elseif (CMAKE_BUILD_TYPE MATCHES Release)
		add_compile_options(-O3 -g -traceback)
	elseif (CMAKE_BUILD_TYPE MATCHES Fast)
		add_compile_options(-O3 -march=native -mtune=native)
	endif()
else()
	add_compile_options(-Wall -Wno-unused -Wno-unused-dummy-argument)
	add_link_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fopenmp>)
	add_compile_options(-g)
	if (CMAKE_BUILD_TYPE MATCHES Debug)
		add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-C>)
		add_compile_options(-O0 -g -ggdb -fbacktrace
			-ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all,no-array-temps)
	elseif (CMAKE_BUILD_TYPE MATCHES Profile)
		add_compile_options(-O2 -p -g -shared-libgcc)
	elseif (CMAKE_BUILD_TYPE MATCHES Release)
		add_compile_options(-O3 -g -fbacktrace)
	elseif (CMAKE_BUILD_TYPE MATCHES Fast)
		add_compile_options(-O3 -march=native -mtune=native)
	endif()
endif()

# RKF45 from https://people.sc.fsu.edu/~jburkardt/
# licenced under LGPL (see SRC/contrib/LICENSE.rkf45)
add_library ( rkf45 SHARED
	SRC/contrib/rkf45.f90
)

add_library ( simple SHARED
	SRC/canonical_coordinates_mod.f90
	SRC/nctools_module.f90
#	SRC/odeint_allroutines.f   # to use this, copy local version into SRC
	SRC/odeint_rkf45.f90       # interface to contrib/rkf45.f90
	SRC/contrib/minpack.f90
	SRC/magfie.f90
	SRC/boozer_converter.f90
	SRC/binsrc.f90
	SRC/plag_coeff.f90
	SRC/chamb_m.f90
	SRC/sub_alpha_lifetime_can.f90
	SRC/vmecinm_m.f90
	SRC/spline_vmec_data.f90
	SRC/spl_three_to_five.f90
	SRC/new_vmec_allocation_stuff.f90
	SRC/get_canonical_coordinates.f90
	SRC/testing.f90
	SRC/field_can.f90
	SRC/orbit_symplectic.f90
	SRC/orbit_symplectic_quasi.f90
	SRC/util.f90
	SRC/simple.f90
	SRC/bench.f90
	SRC/parse_ants.f90
	SRC/zzg.f90
	SRC/rng.f
	SRC/collis_alphas.f90
	SRC/params.f90
)
link_libraries(simple)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	target_link_libraries(simple rkf45 netcdf netcdff)
else()
	target_link_libraries(simple
		rkf45 netcdf netcdff ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
	)
endif()

add_executable (simple.x
	SRC/sorting.f90
	SRC/check_orbit_type.f90
	SRC/find_bminmax.f90
	SRC/main.f90
)

# add_executable (simple_mpi.x
# 	SRC/sorting.f90
# 	SRC/check_orbit_type.f90
# 	SRC/find_bminmax.f90
# 	SRC/main_mpi.f90
# )
# target_link_libraries(simple_mpi.x ${MPI_Fortran_LIBRARIES})

add_executable (test_coord_trans.x SRC/test_coord_trans.f90)
add_executable (test_sympl_tok.x SRC/test_sympl.f90)
add_executable (test_sympl.x SRC/test_sympl_stell.f90)
add_executable (test_magfie.x SRC/test_magfie.f90)
add_executable (test_poincare1.x SRC/test_poincare1.f90)
add_executable (test_poiplot_classification.x
	SRC/test_poiplot_classification.f90
)
add_executable (alpha_lifetime_vmec_c.x
	SRC/binsrc.f90
	SRC/alpha_lifetime_c.f90
)
add_executable (alpha_lifetime_old.x
	SRC/binsrc.f90
	SRC/alpha_lifetime.f90
)
add_executable (test_parse_ants.x SRC/test_parse_ants.f90)
add_executable (test_boozer.x SRC/test_boozer.f90)
add_executable (test_orbits.x SRC/test_orbits_vmec.f90)
add_executable (test_collis.x SRC/test_collis.f90)

add_executable(test_vmec.x EXCLUDE_FROM_ALL SRC/test_vmec.f90)
target_include_directories(test_vmec.x PRIVATE
	$ENV{HOME}/src/STELLOPT/LIBSTELL/Release)
add_library(stell SHARED IMPORTED)
set_target_properties(stell PROPERTIES IMPORTED_LOCATION
	$ENV{HOME}/src/STELLOPT/LIBSTELL/Release/libstell.so)
target_link_libraries(test_vmec.x stell)

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

# Vars
set(f2py_module_name "_pysimple")
set(f2py_module_c "${f2py_module_name}module.c")
set(fortran_src_files
	f90wrap_simple.f90
	f90wrap_orbit_symplectic.f90
	f90wrap_field_can.f90
)
set(generated_module_file "${f2py_module_name}.${Python_SOABI}")
set(f2py_wrappers "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}")

# Generate sources
add_custom_target(
  genpyf
  DEPENDS "${f2py_wrappers}"
)

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
	COMMAND f90wrap -m pysimple
	${CMAKE_CURRENT_LIST_DIR}/SRC/simple.f90
	${CMAKE_CURRENT_LIST_DIR}/SRC/orbit_symplectic.f90
	${CMAKE_CURRENT_LIST_DIR}/SRC/field_can.f90
)

#
# list(APPEND f2py_wrappers "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers2.f90")

# # Set up target
add_library(${f2py_module_name} SHARED
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
