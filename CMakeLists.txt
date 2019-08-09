# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.0.0)
project (NEO-ORB Fortran)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

execute_process(COMMAND "nf-config" --prefix
		OUTPUT_VARIABLE NFPREFIX)

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
	add_compile_options(-mkl -qopenmp -warn all,nounused -stand f95)
	add_link_options(-mkl -qopenmp)
	if (CMAKE_BUILD_TYPE MATCHES Debug)
		add_compile_options(-O0 -g)
	elseif (CMAKE_BUILD_TYPE MATCHES Profile)
		add_compile_options(-O2 -g -shared-intel -debug inline-debug-info
			-D TBB_USE_THREADING_TOOLS -qopenmp-link dynamic -parallel-source-info=2)
	else()
		add_compile_options(-O3 -march=haswell -mtune=haswell)
	endif()
else()
	add_compile_options(-std=f95 -Wall -Wno-unused -Wno-unused-dummy-argument -fopenmp)
	add_link_options(-fopenmp)
	if (CMAKE_BUILD_TYPE MATCHES Debug)
		add_compile_options(-O0 -g -ggdb -C -p -fbacktrace 
			-ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all,no-array-temps)
	elseif (CMAKE_BUILD_TYPE MATCHES Profile)
		add_compile_options(-O2 -g -shared-libgcc)
	else()
		add_compile_options(-O3 -march=ivybridge -mtune=ivybridge)
	endif()
endif()

add_library ( neo_orb SHARED

	SRC/canonical_coordinates_mod.f90
	SRC/nctools_module.f90
	SRC/odeint_allroutines.f
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
	SRC/minpack.f90
	SRC/common.f90
	SRC/neo_orb.f90
)

add_executable (neo-orb.x
#	SRC/canonical_coordinates.f90
	SRC/main.f90
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

add_executable (alpha_lifetime.x
	SRC/zzg.f90
	SRC/binsrc.f90
	SRC/alpha_lifetime.f90
)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	target_link_libraries(neo_orb netcdf netcdff)
else()
	target_link_libraries(neo_orb netcdf netcdff lapack blas)
endif()

target_link_libraries(neo-orb.x neo_orb)
target_link_libraries(test_sympl_tok.x neo_orb)
target_link_libraries(test_sympl.x neo_orb)
target_link_libraries(test_magfie.x neo_orb)
target_link_libraries(test_poincare1.x neo_orb)
target_link_libraries(test_poiplot_classification.x neo_orb)
target_link_libraries(alpha_lifetime_vmec_c.x neo_orb)
target_link_libraries(alpha_lifetime.x neo_orb)
