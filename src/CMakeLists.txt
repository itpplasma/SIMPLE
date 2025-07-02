add_subdirectory(contrib)

set(SOURCES
    coordinates/coordinates.f90
    field/field_base.f90
    field/field_coils.f90
    field/field_vmec.f90
    field/field_gvec.f90
    field.f90
    field/field_can_base.f90
    field/field_can_test.f90
    field/field_can_flux.f90
    field/field_can_boozer.f90
    field/field_can_meiss.f90
    field/psi_transform.f90
    field/field_can_albert.f90
    field_can.f90
    magfie.f90
    boozer_converter.f90
    chamb_m.f90
    sub_alpha_lifetime_can.f90
    vmecinm_m.f90
    spline_vmec_data.f90
    new_vmec_allocation_stuff.f90
    get_canonical_coordinates.f90
    orbit_symplectic_base.f90
    orbit_symplectic_quasi.f90
    orbit_symplectic.f90
    util.f90
    samplers.f90
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


add_library (simple STATIC ${SOURCES})

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	target_link_libraries(simple PRIVATE CONTRIB netcdf netcdff)
else()
	target_link_libraries(simple PUBLIC
		netcdf netcdff CONTRIB BLAS::BLAS LAPACK::LAPACK
	)
endif()

target_link_libraries(simple PUBLIC
    LIBNEO::neo
)
