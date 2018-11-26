# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.11)
project (NEO-ORB Fortran)

include_directories ($ENV{NETCDFF_INCLUDE} /usr/include)
link_directories ($ENV{NETCDF_LIB} $ENV{NETCDFF_LIB})

#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f95 -fall-intrinsics")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -Wall")

add_executable (neo-orb.x
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
	SRC/minpack.f90
	SRC/field_can.f90
	SRC/orbit_symplectic.f90
	SRC/common.f90

	SRC/canonical_coordinates.f90
)

add_executable (orbit-symplectic-test.x
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
 	SRC/minpack.f90
 	SRC/field_can.f90
	SRC/orbit_symplectic.f90
	SRC/common.f90

 	SRC/orbit_symplectic_test.f90
)

add_executable (alpha_lifetime_vmec_c.x
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
	SRC/zzg.f90
	SRC/field_can.f90
	SRC/orbit_symplectic.f90
	SRC/binsrc.f90
	SRC/common.f90

	SRC/alpha_lifetime_c.f90
)


add_executable (test_magfie.x
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
 	SRC/minpack.f90
 	SRC/field_can.f90
	SRC/orbit_symplectic.f90
	SRC/common.f90

 	SRC/test_magfie.f90
)

add_executable (test_poincare1.x
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
	SRC/minpack.f90
	SRC/field_can.f90
	SRC/orbit_symplectic.f90
	SRC/common.f90

	SRC/test_poincare1.f90
)


target_link_libraries(neo-orb.x netcdf netcdff lapack)
target_link_libraries(orbit-symplectic-test.x netcdf netcdff lapack)
target_link_libraries(alpha_lifetime_vmec_c.x netcdf netcdff lapack)
target_link_libraries(test_magfie.x netcdf netcdff lapack)
target_link_libraries(test_poincare1.x netcdf netcdff lapack)
