add_subdirectory(contrib)

set(SOURCES
    coordinates/coordinates.f90
    field/field_base.f90
    field/field_coils.f90
    field/field_vmec.f90
    field.f90
    field/field_can_base.f90
    field/field_can_test.f90
    field/field_can_flux.f90
    field/field_can_boozer.f90
    field/field_can_meiss.f90
    field/psi_transform.f90
    field/field_can_albert.f90
    field_can.f90
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
    orbit_symplectic_base.f90
    orbit_symplectic_quasi.f90
    orbit_symplectic.f90
    util.f90
    cut_detector.f90
    classification.f90
    simple.f90
    parse_ants.f90
    collis_alphas.f90
    callback.f90
    params.f90
    sorting.f90
    check_orbit_type.f90
    find_bminmax.f90
    simple_main.f90
)


add_library (simple SHARED ${SOURCES})

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	target_link_libraries(simple PRIVATE CONTRIB netcdf netcdff)
else()
	target_link_libraries(simple PUBLIC
		netcdf netcdff CONTRIB ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
	)
endif()

### Find HDF5
find_package(HDF5 COMPONENTS Fortran Fortran_HL)
if (NOT ${HDF5_FOUND})
	message(FATAL_ERROR "HDF5 not found, please set environment variable HDF5_ROOT to the correct path")
endif ()

if (APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
  message(STATUS "Apple M1 detected, fixing paths to Fortran HDF5 libraries")
  set(HDF5_LIBRARIES /opt/homebrew/lib/libhdf5_fortran.dylib)
  set(HDF5_HL_LIBRARIES /opt/homebrew/lib/libhdf5_hl_fortran.dylib)
endif()

# Set variables in parent scope, so find_package is not needed there
set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS} PARENT_SCOPE)
set(HDF5_LIBRARIES ${HDF5_LIBRARIES} PARENT_SCOPE)
set(HDF5_HL_LIBRARIES ${HDF5_HL_LIBRARIES} PARENT_SCOPE)


message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")
message(STATUS "HDF5 HL libraries: ${HDF5_HL_LIBRARIES}")

include_directories(${HDF5_INCLUDE_DIRS})

add_library(hdf5_tools
    libneo/src/hdf5_tools/hdf5_tools.f90
	libneo/src/hdf5_tools/hdf5_tools_f2003.f90
)
target_link_libraries(simple PUBLIC $<TARGET_OBJECTS:hdf5_tools>)
target_link_libraries(simple PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

add_subdirectory(libneo/src/interpolate interpolate)
target_link_libraries(simple PUBLIC $<TARGET_OBJECTS:interpolate>)

add_subdirectory(libneo/src/polylag neo_polylag)
target_link_libraries(simple PUBLIC $<TARGET_OBJECTS:neo_polylag>)

add_subdirectory(libneo/src/hdf5_tools)
target_link_libraries(simple PUBLIC hdf5_tools)

add_subdirectory(libneo/src/interpolate interpolate)
target_link_libraries(simple PUBLIC $<TARGET_OBJECTS:interpolate>)

add_subdirectory(libneo/src/field neo_field)
add_dependencies(neo_field neo_polylag interpolate hdf5_tools)
target_link_libraries(simple PUBLIC $<TARGET_OBJECTS:neo_field>)