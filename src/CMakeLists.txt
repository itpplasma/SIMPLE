set(SOURCES
    lapack_interfaces.f90
    sorting_mod.f90
    coordinates/coordinates.f90
    coordinates/stencil_utils.f90
    coordinates/array_utils.f90
    field/field_base.f90
    field/field_analytical_gs.f90
    field/field_coils.f90
    field/field_vmec.f90
    field/field_geoflux.f90
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
    magfie_can_boozer.f90
    boozer_converter.F90
    chamb_m.f90
    sub_alpha_lifetime_can.f90
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
    netcdf_orbit_output.f90
    callback.f90
    params.f90
    sorting.f90
    check_orbit_type.f90
    find_bminmax.f90
    timing.f90
    simple_main.f90
)

# Add diagnostic sources
list(APPEND SOURCES
    diag/diag_meiss.f90
    diag/diag_albert.f90
    diag/diag_newton.f90
    diag/diag_orbit.f90
)

# Add GVEC-specific sources if enabled
if(ENABLE_GVEC)
    list(APPEND SOURCES
        field/field_gvec.f90
        field/vmec_field_adapter.f90
    )
endif()

add_library (simple STATIC ${SOURCES})

set(_magfie_stamp "${CMAKE_BINARY_DIR}/magfie_modules.stamp")
add_custom_command(
    OUTPUT ${_magfie_stamp}
    COMMAND ${CMAKE_COMMAND} -E touch ${_magfie_stamp}
    DEPENDS magfie
)
add_custom_target(prepare_magfie_modules DEPENDS ${_magfie_stamp})

set_source_files_properties(field/field_geoflux.f90
    PROPERTIES OBJECT_DEPENDS ${_magfie_stamp})
set_source_files_properties(coordinates/coordinates.f90
    PROPERTIES OBJECT_DEPENDS ${_magfie_stamp})


# Link pyplot library to SIMPLE
target_link_libraries(simple PUBLIC pyplot)

# Apply SIMPLE-specific compile options
target_compile_options(simple PRIVATE ${SIMPLE_COMPILE_OPTIONS})

# Apply trampoline error flags only to SIMPLE library (not subprojects like fortplot)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(simple PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wtrampolines>
        $<$<COMPILE_LANGUAGE:Fortran>:-Werror=trampolines>
    )
endif()

target_link_libraries(simple PUBLIC
	netcdf netcdff CONTRIB BLAS::BLAS LAPACK::LAPACK
)

target_link_libraries(simple PUBLIC
    LIBNEO::neo
    LIBNEO::magfie
)

add_dependencies(simple neo prepare_magfie_modules)

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

# Link fortplot
target_link_libraries(simple PUBLIC
    fortplot::fortplot
)

# Link coverage flags to library when coverage is enabled
if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE MATCHES "Debug|Profile")
    target_link_libraries(simple PUBLIC coverage_flags)
endif()
