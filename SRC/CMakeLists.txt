add_subdirectory(contrib)

set(SOURCES
    canonical_coordinates_mod.f90
    nctools_module.f90
    #    odeint_allroutines.f   # to use this, copy local version into SRC
    odeint_rkf45.f90       # interface to contrib/rkf45.f90
    magfie.f90
    boozer_converter.f90
    binsrc.f90
    coordinates/coordinates.f90
    plag_coeff.f90
    chamb_m.f90
    sub_alpha_lifetime_can.f90
    vmecinm_m.f90
    spline_vmec_data.f90
    spl_three_to_five.f90
    new_vmec_allocation_stuff.f90
    get_canonical_coordinates.f90
    interpolate.f90
    magfie/magfie_base.f90
    magfie/magfie_coils.f90
    magfie/magfie_vmec.f90
    magfie/magfie.f90
    field/field_vmec.f90
    field/field_can_base.f90
    field/field_can_test.f90
    field/field_can_flux.f90
    field/field_can_boozer.f90
    field/field_can_meiss.f90
    field/psi_transform.f90
    field/field_can_albert.f90
    field_can.f90
    orbit_symplectic_base.f90
    orbit_symplectic_quasi.f90
    orbit_symplectic.f90
    util.f90
    cut_detector.f90
    classification.f90
    simple.f90
    parse_ants.f90
    collis_alphas.f90
    params.f90
    sorting.f90
    check_orbit_type.f90
    find_bminmax.f90
    callback.f90
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

add_subdirectory(${CMAKE_SOURCE_DIR}/../libneo/src/polylag neo_polylag)
target_link_libraries(simple PUBLIC $<TARGET_OBJECTS:neo_polylag>)

add_subdirectory(${CMAKE_SOURCE_DIR}/../libneo/src/field neo_field)
add_dependencies(neo_field neo_polylag)
target_link_libraries(simple PUBLIC $<TARGET_OBJECTS:neo_field>)
