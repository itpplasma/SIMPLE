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
        coordinates/scaled_chartmap_coordinates.f90
    coordinates/cartesian_coordinates.f90
    field/field_base.f90
    field/field_coils.f90
    field/field_vmec.f90
    field/field_splined.f90
    field/boozer_chartmap_io.f90
    field/field_boozer_chartmap.f90
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
    orbit_symplectic_euler1.f90
    orbit_rk_core.f90
    orbit_symplectic.f90
    orbit_cpp.f90
    field_pauli_cart.f90
    pauli_gc_drift.f90
    orbit_cpp_pauli.f90
    orbit_cpp_vmec_metric.f90
    orbit_cpp_canonical.f90
    orbit_full_provider.f90
    orbit_full_mock_cart.f90
    orbit_full_mock_cyl.f90
    orbit_full.f90
    orbit_full_device.f90
    util.F90
    samplers.f90
    cut_detector.f90
    classification.f90
    simple.f90
    parse_ants.f90
    profiles.f90
    collis_alphas.f90
    netcdf_orbit_output.f90
    netcdf_results_output.f90
    callback.f90
    params.f90
    params_wrapper.f90
    sorting.f90
    check_orbit_type.f90
    find_bminmax.f90
    diag_counters.f90
    restart.f90
    progress_monitor.f90
    timing.f90
    simple_gpu.f90
    simple_main.f90
    simple_entry.f90
    )

    if(SIMPLE_ENABLE_CGAL)
        add_library(stl_wall_cgal STATIC wall/stl_wall_cgal.cpp)
        target_include_directories(stl_wall_cgal PUBLIC ${CGAL_INCLUDE_DIR})

        # CGAL_DISABLE_GMP makes CGAL's bundled CORE use Boost.Multiprecision
        # cpp_int instead of GMP-backed mpz_int, so the wall feature needs only
        # Boost headers. Without it, CORE pulls in GMP (CORE::BigInt is
        # boost::multiprecision::mpz_int) and the build fails wherever gmp.h is
        # absent, which is most HPC nodes (only libgmp.so, no -dev package).
        target_compile_definitions(stl_wall_cgal PUBLIC CGAL_DISABLE_GMP)
        target_link_libraries(stl_wall_cgal PUBLIC Boost::headers)
    endif()

# Add diagnostic sources
list(APPEND SOURCES
    diag/diag_meiss.f90
    diag/diag_albert.f90
    diag/diag_newton.f90
    diag/diag_orbit.f90
)

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
    netcdf netcdff BLAS::BLAS LAPACK::LAPACK
)

target_link_libraries(simple PUBLIC
    LIBNEO::neo
)

target_link_libraries(simple PUBLIC
    fortnum
)

# For NVHPC: add built deps includes first to avoid gfortran-compiled modules in /usr/include
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_include_directories(simple BEFORE PUBLIC
        ${CMAKE_BINARY_DIR}/deps/include
    )
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
