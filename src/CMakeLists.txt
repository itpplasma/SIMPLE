    set(SOURCES
        ${CMAKE_BINARY_DIR}/version.f90
        lapack_interfaces.f90
        sorting_mod.f90
        wall/chartmap_metadata.f90
        wall/stl_wall_intersection.F90
        coordinates/coordinates.f90
        coordinates/stencil_utils.f90
        coordinates/array_utils.f90
        coordinates/reference_coordinates.f90
        coordinates/coordinate_scaling.f90
    coordinates/cartesian_coordinates.f90
    field/field_base.f90
    field/field_coils.f90
    field/field_vmec.f90
    field/field_splined.f90
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
    magfie_wrapper.f90
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
    netcdf_results_output.f90
    callback.f90
    params.f90
    params_wrapper.f90
    sorting.f90
    check_orbit_type.f90
    find_bminmax.f90
        timing.f90
        simple_main.f90
        simple_entry.f90
    )

    if(SIMPLE_ENABLE_CGAL)
        add_library(stl_wall_cgal STATIC wall/stl_wall_cgal.cpp)
        # Header-only usage with Simple_cartesian kernel (no GMP/MPFR needed)
        if(TARGET CGAL::CGAL)
            target_link_libraries(stl_wall_cgal PUBLIC CGAL::CGAL)
        else()
            # Fetched CGAL: just need include path (header-only)
            target_include_directories(stl_wall_cgal PUBLIC ${CGAL_INCLUDE_DIR})
        endif()
    endif()

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

    if(SIMPLE_ENABLE_CGAL)
        target_link_libraries(simple PUBLIC stl_wall_cgal)
    endif()

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
)

# For NVHPC: add built deps includes first to avoid gfortran-compiled modules in /usr/include
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_include_directories(simple BEFORE PUBLIC
        ${CMAKE_BINARY_DIR}/deps/include
    )
endif()

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

# Link fortplot (when available)
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_link_libraries(simple PUBLIC
        fortplot::fortplot
    )
endif()

# Link coverage flags to library when coverage is enabled
if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE MATCHES "Debug|Profile")
    target_link_libraries(simple PUBLIC coverage_flags)
endif()
