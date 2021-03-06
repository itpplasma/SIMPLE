# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.12)
project (SIMPLE Fortran)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_program(NF_CONFIG "nf-config")

if (NF_CONFIG)
execute_process(COMMAND "nf-config" --prefix
		OUTPUT_VARIABLE NFPREFIX)
else()
message(SEND_ERROR "nf-config not found. Please install libnetcdff-dev")
endif()

string(STRIP ${NFPREFIX} NFPREFIX)
set(NFINC ${NFPREFIX}/include)
set(NFLIBS ${NFPREFIX}/lib)

message(STATUS "CMake build type: " ${CMAKE_BUILD_TYPE})
message(STATUS "NetCDF include path: " ${NFINC})
message(STATUS "NetCDF lib path: " ${NFLIBS})

include_directories ($ENV{NETCDFF_INCLUDE} ${NFINC})
link_directories ($ENV{NETCDF_LIB} $ENV{NETCDFF_LIB} ${NFLIBS})

add_compile_options(-cpp)
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	add_compile_options(-mkl -qopenmp -warn all,nounused -stand f2008)
	link_libraries("-mkl -qopenmp")
	if (CMAKE_BUILD_TYPE MATCHES Debug)
		add_compile_options(-O0 -g)
	elseif (CMAKE_BUILD_TYPE MATCHES Profile)
		add_compile_options(-O2 -g -shared-intel -debug inline-debug-info
			-D TBB_USE_THREADING_TOOLS -qopenmp-link dynamic -parallel-source-info=2)
	elseif (CMAKE_BUILD_TYPE MATCHES Fast)
		add_compile_options(-O3 -march=skylake -mtune=skylake)
	endif()
else()
	add_compile_options(-std=f2008 -Wall -Wno-unused -Wno-unused-dummy-argument -fopenmp)
	link_libraries("-fopenmp")
	if (CMAKE_BUILD_TYPE MATCHES Debug)
		add_compile_options(-O0 -g -ggdb -C -p -fbacktrace
			-ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all,no-array-temps)
	elseif (CMAKE_BUILD_TYPE MATCHES Profile)
		add_compile_options(-O2 -g -shared-libgcc)
	elseif (CMAKE_BUILD_TYPE MATCHES Fast)
		add_compile_options(-O3 -march=ivybridge -mtune=ivybridge)
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
	SRC/common.f90
	SRC/simple.f90
	SRC/bench.f90
	SRC/parse_ants.f90
)

add_executable (simple.x
	SRC/zzg.f90
	SRC/binsrc.f90
	SRC/main.f90
)

add_executable (test_coord_trans.x
	SRC/zzg.f90
	SRC/binsrc.f90
	SRC/test_coord_trans.f90
)

add_executable (test_sympl_tok.x
  	SRC/test_sympl.f90
)

add_executable (test_sympl.x
  	SRC/test_sympl_stell.f90
)

add_executable (test_magfie.x
  	SRC/test_magfie.f90
)

add_executable (test_poincare1.x
	SRC/test_poincare1.f90
)

add_executable (test_poiplot_classification.x
	SRC/test_poiplot_classification.f90
)

add_executable (alpha_lifetime_vmec_c.x
	SRC/zzg.f90
	SRC/binsrc.f90
	SRC/alpha_lifetime_c.f90
)

add_executable (alpha_lifetime_old.x
	SRC/zzg.f90
	SRC/binsrc.f90
	SRC/alpha_lifetime.f90
)

add_executable (test_parse_ants.x
	SRC/parse_ants.f90
	SRC/test_parse_ants.f90
)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	target_link_libraries(simple rkf45 netcdf netcdff)
else()
	target_link_libraries(simple rkf45 netcdf netcdff lapack blas)
endif()

target_link_libraries(simple.x simple)
target_link_libraries(test_coord_trans.x simple)
target_link_libraries(test_sympl_tok.x simple)
target_link_libraries(test_sympl.x simple)
target_link_libraries(test_magfie.x simple)
target_link_libraries(test_poincare1.x simple)
target_link_libraries(test_poiplot_classification.x simple)
target_link_libraries(alpha_lifetime_vmec_c.x simple)
target_link_libraries(alpha_lifetime_old.x simple)
