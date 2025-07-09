add_subdirectory(contrib)

set(SOURCES
    lapack_interfaces.f90
    sorting_mod.f90
    coordinates/coordinates.f90
    coordinates/stencil_utils.f90
    coordinates/array_utils.f90
    field/field_base.f90
    field/field_coils.f90
    field/field_vmec.f90
    field/vmec_field_eval.f90
    field/field_newton.F90
    field.F90
    field/field_can_base.f90
    field/field_can_test.f90
    field/field_can_flux.f90
    field/field_can_boozer.f90
    field/field_can_meiss.f90
    field/psi_transform.f90
    field/field_can_albert.f90
    field_can.f90
    magfie.f90
    boozer_converter.F90
    chamb_m.f90
    sub_alpha_lifetime_can.f90
    spline_vmec_data.f90
    new_vmec_allocation_stuff.f90
    get_canonical_coordinates.F90
    orbit_symplectic_base.f90
    orbit_symplectic_quasi.f90
    orbit_symplectic.f90
    util.F90
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

# Add GVEC-specific sources if enabled
if(ENABLE_GVEC)
    list(APPEND SOURCES
        field/field_gvec.f90
        field/vmec_field_adapter.f90
    )
endif()

add_library (simple STATIC ${SOURCES})

target_link_libraries(simple PUBLIC
	netcdf netcdff CONTRIB BLAS::BLAS LAPACK::LAPACK
)

target_link_libraries(simple PUBLIC
    LIBNEO::neo
)

# Conditionally link GVEC if enabled
if(ENABLE_GVEC)
    target_include_directories(simple PRIVATE
        ${gvec_BINARY_DIR}/include
    )
    target_link_libraries(simple PUBLIC
        gveclib
    )
    target_compile_definitions(simple PUBLIC GVEC_AVAILABLE)
endif()
