# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.12)
project (SIMPLE Fortran C)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
#    set(CMAKE_BUILD_TYPE Debug)
endif()

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
add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008>)
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
		add_compile_options(-O0 -g -ggdb -C -fbacktrace
			-ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all,no-array-temps)
	elseif (CMAKE_BUILD_TYPE MATCHES Profile)
		add_compile_options(-O2 -p -g -shared-libgcc)
	elseif (CMAKE_BUILD_TYPE MATCHES Release)
		add_compile_options(-O3 -g -fbacktrace)
	elseif (CMAKE_BUILD_TYPE MATCHES Fast)
		add_compile_options(-O3 -march=native -mtune=native)
	endif()
endif()

if (CMAKE_C_COMPILER_ID MATCHES Clang)
	set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
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

add_executable (simple_mpi.x
	SRC/sorting.f90
	SRC/check_orbit_type.f90
	SRC/find_bminmax.f90
	SRC/main_mpi.f90
)
target_link_libraries(simple_mpi.x ${MPI_Fortran_LIBRARIES})

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
find_package(Python3 REQUIRED
  COMPONENTS Interpreter Development NumPy)

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
  COMMAND "${Python3_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

include_directories(
    BEFORE
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
    ${F2PY_INCLUDE_DIR}
)

message(STATUS ${Python3_INCLUDE_DIRS})
message(STATUS ${F2PY_INCLUDE_DIR})
message(STATUS ${Python3_NumPy_INCLUDE_DIRS})

# Vars
set(f2py_module_name "pysimple")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/SRC/pysimple.f90")
set(fortran_pyf_file "${CMAKE_SOURCE_DIR}/SRC/pysimple.pyf")
set(f2py_module_c "${f2py_module_name}module.c")
set(generated_module_file "${f2py_module_name}${Python3_SOABI}")

# Generate sources
add_custom_target(
  genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
  )

set(f2py_wrappers "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}module.c")
list(APPEND f2py_wrappers "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers2.f90")

add_custom_command(
  OUTPUT ${f2py_wrappers}
  COMMAND ${Python3_EXECUTABLE}  -m "numpy.f2py"
                   "${fortran_pyf_file}"
                   -m "${f2py_module_name}"
                   --lower # Important
  DEPENDS "${fortran_src_file}" # Fortran source
  )

# Set up target
add_library(${f2py_module_name} SHARED
	${f2py_wrappers}
  "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
  "${fortran_src_file}" # Fortran source(s)
  )

# Depend on sources
add_dependencies(${f2py_module_name} genpyf)

set_target_properties(
    ${f2py_module_name}
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "${f2py_module_name}"
		SUFFIX ".so"
        LINKER_LANGUAGE C
    )
