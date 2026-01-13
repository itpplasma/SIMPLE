# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.22)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project (SIMPLE)
enable_language(C Fortran)

# Optional: enable CGAL-based STL wall intersection support.
option(SIMPLE_ENABLE_CGAL "Enable CGAL STL wall intersection support" OFF)

# Get version from git describe (format: v1.5.2-42-g62fed6e[-dirty])
execute_process(
    COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE SIMPLE_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT SIMPLE_VERSION)
    set(SIMPLE_VERSION "unknown")
endif()
message(STATUS "SIMPLE version: ${SIMPLE_VERSION}")

# Generate version.f90 from template
configure_file(
    ${CMAKE_SOURCE_DIR}/src/version.f90.in
    ${CMAKE_BINARY_DIR}/version.f90
    @ONLY
)

add_compile_options(-g)

# Compiler-specific flags
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fbacktrace>)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-traceback>)
    # Keep denormals and IEEE division semantics for numerical robustness.
    # Needed across SIMPLE and Fortran dependencies (e.g., libneo ODE integrators).
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Mnoflushz>)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Kieee>)
endif()

# Disable executable stack for security (trampolines/closures create execstack)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # Apple's linker doesn't support -z,noexecstack (uses different security model)
    if(NOT APPLE)
        add_link_options(-Wl,-z,noexecstack)
    endif()
endif()

# Note: Trampoline warnings are applied per-target to SIMPLE code only,
# not globally, to avoid affecting subprojects like fortplot

include(Util)

set	(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
	Release Debug Profile Fast
)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE
	STRING "Choose the type of build." FORCE)
endif()

option(ENABLE_PYTHON_INTERFACE "Enables the Python-Wrapper." ON)
option(SIMPLE_TESTING "Enable testing." ON)
option(ENABLE_OPENMP "Enable OpenMP compiler flags." ON)
option(ENABLE_GVEC "Enable GVEC field support (experimental)" OFF)
option(ENABLE_COVERAGE "Enable code coverage analysis (Debug/Profile builds only)" OFF)
option(SIMPLE_DETERMINISTIC_FP "Disable fast-math for reproducible floating-point" OFF)
option(SIMPLE_ENABLE_PYTHON_TOOLS "Enable Python helpers (tests/data generation)" ON)
option(SIMPLE_ENABLE_OPENACC "Enable OpenACC offload (NVHPC only)" OFF)
option(SIMPLE_ENABLE_DEBUG_OUTPUT "Enable debug output in hot loops" OFF)

if(SIMPLE_ENABLE_DEBUG_OUTPUT)
    add_compile_definitions(SIMPLE_ENABLE_DEBUG_OUTPUT)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    # NVHPC/nvfortran: disable all Python integration/tools by default.
    if(ENABLE_PYTHON_INTERFACE)
        message(STATUS "Python interface disabled for NVHPC compiler")
    endif()
    set(ENABLE_PYTHON_INTERFACE OFF CACHE BOOL "Enables the Python-Wrapper." FORCE)
    if(SIMPLE_ENABLE_PYTHON_TOOLS)
        message(STATUS "Python tools disabled for NVHPC compiler")
    endif()
    set(SIMPLE_ENABLE_PYTHON_TOOLS OFF CACHE BOOL
        "Enable Python helpers (tests/data generation)" FORCE)
    add_compile_definitions(SIMPLE_NVHPC)

    if(SIMPLE_ENABLE_OPENACC)
        message(STATUS "OpenACC enabled for NVHPC compiler")
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-acc>
            $<$<COMPILE_LANGUAGE:Fortran>:-gpu=ccnative>
            $<$<COMPILE_LANGUAGE:C>:-acc>
            $<$<COMPILE_LANGUAGE:C>:-gpu=ccnative>
        )
        add_link_options(
            -acc
            -gpu=ccnative
        )
    endif()
endif()

# Conditionally fetch GVEC if enabled
if(ENABLE_GVEC)
    message(STATUS "GVEC support enabled - fetching GVEC library")
    include(FetchContent)
    FetchContent_Declare(
      gvec
      GIT_REPOSITORY https://gitlab.mpcdf.mpg.de/calbert/gvec.git
      GIT_TAG        develop
    )
    set(GVEC_FOR_SIMPLE TRUE CACHE BOOL "Build minimal GVEC for SIMPLE")
    FetchContent_MakeAvailable(gvec)
else()
    message(STATUS "GVEC support disabled")
endif()

set(CMAKE_MACOSX_RPATH 1)

find_program(NC_CONFIG "nf-config")

if (NC_CONFIG)
execute_process(COMMAND nf-config --includedir
                OUTPUT_VARIABLE NETCDFINCLUDE_DIR)
execute_process(COMMAND nc-config --libdir
				OUTPUT_VARIABLE NETCDFLIB_DIR)
execute_process(COMMAND nf-config --flibs
                OUTPUT_VARIABLE NETCDF_FLIBS)
else()
message(SEND_ERROR "nc-config not found. Please install libnetcdff-dev")
endif()

string(STRIP ${NETCDFINCLUDE_DIR} NETCDFINCLUDE_DIR)
string(STRIP ${NETCDFLIB_DIR} NETCDFLIB_DIR)
string(STRIP ${NETCDF_FLIBS} NETCDF_FLIBS)

message(STATUS "NetCDF include path: " ${NETCDFINCLUDE_DIR})
message(STATUS "NetCDF lib path: " ${NETCDFLIB_DIR})
message(STATUS "NetCDF Fortran libs: " ${NETCDF_FLIBS})
string(REPLACE " " ";" NETCDF_FLIBS "${NETCDF_FLIBS}")

# For NVHPC compiler: avoid global include of /usr/include as it contains
# gfortran-compiled HDF5 modules incompatible with nvfortran
# NetCDF includes will be added target-specifically later
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    include_directories(${NETCDFINCLUDE_DIR})
endif()
link_directories(${NETCDFLIB_DIR})
add_link_options(${NETCDF_FLIBS})

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    # Use NVHPC-provided BLAS/LAPACK to avoid pulling in libgomp via system OpenBLAS.
    get_filename_component(_nvhpc_bin_dir "${CMAKE_Fortran_COMPILER}" DIRECTORY)
    get_filename_component(_nvhpc_compilers_dir "${_nvhpc_bin_dir}/.." REALPATH)
    set(_nvhpc_lib_dir "${_nvhpc_compilers_dir}/lib")

    find_library(NVHPC_BLAS_LIBRARY
        NAMES blas_lp64
        PATHS "${_nvhpc_lib_dir}"
        NO_DEFAULT_PATH
        REQUIRED
    )
    find_library(NVHPC_LAPACK_LIBRARY
        NAMES lapack_lp64
        PATHS "${_nvhpc_lib_dir}"
        NO_DEFAULT_PATH
        REQUIRED
    )

    if(NOT TARGET BLAS::BLAS)
        add_library(BLAS::BLAS INTERFACE IMPORTED)
    endif()
    set_target_properties(BLAS::BLAS PROPERTIES
        INTERFACE_LINK_LIBRARIES "${NVHPC_BLAS_LIBRARY}"
    )

    if(NOT TARGET LAPACK::LAPACK)
        add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    endif()
    set_target_properties(LAPACK::LAPACK PROPERTIES
        INTERFACE_LINK_LIBRARIES "${NVHPC_LAPACK_LIBRARY}"
    )
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()

find_or_fetch(libneo)

# CGAL (optional, for STL wall intersections)
# Uses Simple_cartesian<double> kernel - header-only, no GMP/MPFR needed
if(SIMPLE_ENABLE_CGAL)
    find_package(CGAL QUIET)
    if(CGAL_FOUND)
        message(STATUS "CGAL enabled (system): ${CGAL_VERSION}")
        set(CGAL_INCLUDE_DIR "")  # Use CGAL::CGAL target
    else()
        message(STATUS "CGAL not found on system, fetching header-only version...")
        include(FetchContent)
        FetchContent_Declare(
            cgal
            URL https://github.com/CGAL/cgal/releases/download/v5.6.2/CGAL-5.6.2-library.tar.xz
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        )
        FetchContent_Populate(cgal)
        set(CGAL_INCLUDE_DIR "${cgal_SOURCE_DIR}/include")
        message(STATUS "CGAL enabled (fetched): 5.6.2 at ${CGAL_INCLUDE_DIR}")
    endif()
    add_compile_definitions(SIMPLE_ENABLE_CGAL)
else()
    message(STATUS "CGAL disabled (SIMPLE_ENABLE_CGAL=OFF)")
endif()

# Fetch fortplot for diagnostic plotting (disabled for nvfortran due to compatibility issues)
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    include(FetchContent)
    FetchContent_Declare(
        fortplot
        GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
        GIT_TAG feabb592faf39bac8167b9193413912f2bb8ca48
    )
    FetchContent_MakeAvailable(fortplot)
    add_compile_definitions(USE_FORTPLOT)
else()
    message(STATUS "Fortplot disabled for NVHPC compiler (compatibility issues)")
endif()

if (SIMPLE_TESTING)
    message(STATUS "Unit Tests enabled!")
endif()

message(STATUS "CMake build type: " ${CMAKE_BUILD_TYPE})

include_directories ($ENV{NETCDFF_INCLUDE} ${NFINC})
link_directories ($ENV{NETCDF_LIB} $ENV{NETCDFF_LIB} ${NFLIBS})

# Set SIMPLE-specific compiler flags after GVEC to avoid affecting it
# Compiler-specific flags (not applied to subprojects)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(SIMPLE_COMPILE_OPTIONS
        -Wall -Wextra -fPIC -g
        $<$<COMPILE_LANGUAGE:Fortran>:-Wimplicit-interface>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-external-argument-mismatch>
        $<$<COMPILE_LANGUAGE:Fortran>:-fmax-errors=5>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-132>
        $<$<COMPILE_LANGUAGE:Fortran>:-Werror=line-truncation>
        $<$<COMPILE_LANGUAGE:Fortran>:-fbacktrace>
    )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(SIMPLE_COMPILE_OPTIONS
        -fPIC -g
        $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>
        # Avoid flush-to-zero altering numerics (notably chartmap/Meiss ODE integration).
        $<$<COMPILE_LANGUAGE:Fortran>:-Mnoflushz>
        # Use IEEE division semantics for numerical robustness.
        $<$<COMPILE_LANGUAGE:Fortran>:-Kieee>
        $<$<COMPILE_LANGUAGE:Fortran>:-traceback>
    )
else()
    set(SIMPLE_COMPILE_OPTIONS -fPIC -g)
endif()

if (ENABLE_OPENMP)
	# add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fopenmp>)
	# add_link_options(-fopenmp)
    find_package(OpenMP REQUIRED)
    link_libraries(OpenMP::OpenMP_Fortran)
endif()

# Coverage analysis configuration (only for Debug/Profile builds)
if (ENABLE_COVERAGE)
    if (CMAKE_BUILD_TYPE MATCHES "Debug|Profile")
        message(STATUS "Coverage analysis enabled for ${CMAKE_BUILD_TYPE} build")
        
        # Create interface library for coverage flags to scope them to specific targets
        add_library(coverage_flags INTERFACE)
        target_compile_options(coverage_flags INTERFACE --coverage -fprofile-arcs -ftest-coverage)
        target_link_options(coverage_flags INTERFACE --coverage -lgcov)
        
        # Add custom target for coverage data generation
        find_program(LCOV_PATH lcov)
        if(LCOV_PATH)
            add_custom_target(coverage
                COMMAND ${CMAKE_COMMAND} -E echo "Generating coverage data..."
                COMMAND ${LCOV_PATH} --capture --directory ${CMAKE_BINARY_DIR} --output-file coverage.info
                    --rc branch_coverage=1 
                    --ignore-errors inconsistent 
                    --ignore-errors mismatch 
                    --ignore-errors unused
                COMMAND ${LCOV_PATH} --remove coverage.info
                    "/home/ert/code/libneo/*"
                    "*/libneo/*"
                    "*/test/*"
                    "/usr/*"
                    --output-file coverage_filtered.info
                    --rc branch_coverage=1
                    --ignore-errors mismatch
                    --ignore-errors unused
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                COMMENT "Generating code coverage report"
                VERBATIM
            )
        else()
            message(WARNING "lcov not found - coverage target will not be available")
        endif()
    else()
        message(WARNING "Coverage analysis requires Debug or Profile build type (current: ${CMAKE_BUILD_TYPE})")
    endif()
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
	add_compile_options(-O0 -ggdb -fno-omit-frame-pointer -fno-inline)
	add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-C>)
	add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ffpe-trap=invalid,zero,overflow>)
	add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fbounds-check>)
	add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fcheck=all,no-array-temps>)
	add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-finit-real=nan>)
elseif (CMAKE_BUILD_TYPE MATCHES Profile)
	add_compile_options(-Og -pg)
	add_link_options(-pg)
elseif (CMAKE_BUILD_TYPE MATCHES "Release|Fast")
	# Common optimization flags for Release and Fast builds
	if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
		if(SIMPLE_DETERMINISTIC_FP)
			# Deterministic floating-point: no fast-math, strict contract
			list(APPEND SIMPLE_COMPILE_OPTIONS -O3 -ffp-contract=off -funroll-loops)
			message(STATUS "Deterministic FP mode: disabled fast-math for reproducibility")
		else()
			list(APPEND SIMPLE_COMPILE_OPTIONS -O3 -ffast-math -ffp-contract=fast -funroll-loops)
		endif()

		# Architecture-specific flags for GNU
		if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
			# ARM processors (e.g., Apple Silicon): gfortran doesn't support -march=native
			# Use -mcpu=native instead for ARM
			list(APPEND SIMPLE_COMPILE_OPTIONS -mcpu=native)
		elseif (CMAKE_BUILD_TYPE MATCHES Fast)
			# Fast mode on x86: use native architecture
			list(APPEND SIMPLE_COMPILE_OPTIONS -march=native -mtune=native)
		else()
			# Release on x86-64: use v2 baseline for broad compatibility
			list(APPEND SIMPLE_COMPILE_OPTIONS -march=x86-64-v2 -mtune=native)
		endif()
	elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
		# NVHPC optimization flags
		if(SIMPLE_DETERMINISTIC_FP)
			# Deterministic floating-point: no fast-math, strict IEEE
			list(APPEND SIMPLE_COMPILE_OPTIONS -O3)
			message(STATUS "Deterministic FP mode: disabled fast-math for reproducibility (NVHPC)")
		else()
			list(APPEND SIMPLE_COMPILE_OPTIONS -fast -O3 -Mfprelaxed)
		endif()
	else()
		list(APPEND SIMPLE_COMPILE_OPTIONS -O3)
	endif()
elseif (CMAKE_BUILD_TYPE MATCHES Debug)
	if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
		list(APPEND SIMPLE_COMPILE_OPTIONS -O0 -ggdb -fno-omit-frame-pointer -fno-inline)
		list(APPEND SIMPLE_COMPILE_OPTIONS $<$<COMPILE_LANGUAGE:Fortran>:-C>)
		list(APPEND SIMPLE_COMPILE_OPTIONS $<$<COMPILE_LANGUAGE:Fortran>:-ffpe-trap=invalid,zero,overflow>)
		list(APPEND SIMPLE_COMPILE_OPTIONS $<$<COMPILE_LANGUAGE:Fortran>:-fbounds-check>)
		list(APPEND SIMPLE_COMPILE_OPTIONS $<$<COMPILE_LANGUAGE:Fortran>:-fcheck=all,no-array-temps>)
		list(APPEND SIMPLE_COMPILE_OPTIONS $<$<COMPILE_LANGUAGE:Fortran>:-finit-real=nan>)
	elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
		list(APPEND SIMPLE_COMPILE_OPTIONS -O0 -g -Mbounds -Mchkptr)
	else()
		list(APPEND SIMPLE_COMPILE_OPTIONS -O0 -g)
	endif()
elseif (CMAKE_BUILD_TYPE MATCHES Profile)
	if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
		list(APPEND SIMPLE_COMPILE_OPTIONS -O2 -p)
	elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
		list(APPEND SIMPLE_COMPILE_OPTIONS -O2 -Mprof=ccff)
	else()
		list(APPEND SIMPLE_COMPILE_OPTIONS -O2)
	endif()
else()
    message(FATAL_ERROR "Unknown build type: ${CMAKE_BUILD_TYPE}")
endif()

# Trampoline error flags will be applied per-target to avoid affecting subprojects

# Auto-download pyplot-fortran single file during build
set(PYPLOT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)
set(PYPLOT_MODULE_FILE ${PYPLOT_SOURCE_DIR}/pyplot_module.F90)

# Download pyplot_module.F90 if it doesn't exist
if(NOT EXISTS ${PYPLOT_MODULE_FILE})
    message(STATUS "Downloading pyplot_module.F90...")
    file(MAKE_DIRECTORY ${PYPLOT_SOURCE_DIR})
    file(DOWNLOAD
        "https://raw.githubusercontent.com/jacobwilliams/pyplot-fortran/refs/heads/master/src/pyplot_module.F90"
        ${PYPLOT_MODULE_FILE}
        SHOW_PROGRESS
    )
endif()

# Create pyplot library from downloaded source
add_library(pyplot STATIC ${PYPLOT_MODULE_FILE})
set_target_properties(pyplot PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
target_include_directories(pyplot PUBLIC ${CMAKE_BINARY_DIR}/include)

add_subdirectory(src)
include_directories("${CMAKE_BINARY_DIR}/src")

add_executable (simple.x
	app/simple.f90
)

# Add diagnostic executables
add_executable(diag_meiss.x
    app/simple_diag_meiss.f90
)

add_executable(diag_albert.x
    app/simple_diag_albert.f90
)

add_executable(diag_newton.x
    app/simple_diag_newton.f90
)

add_executable(diag_orbit.x
    app/simple_diag_orbit.f90
)

# Apply SIMPLE-specific compile options to executable
target_compile_options(simple.x PRIVATE ${SIMPLE_COMPILE_OPTIONS})

# Apply trampoline error flags only to SIMPLE executable (not subprojects like fortplot)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(simple.x PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wtrampolines>
        $<$<COMPILE_LANGUAGE:Fortran>:-Werror=trampolines>
    )
endif()

target_link_libraries(simple.x PRIVATE simple)
set_target_properties(simple.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Configure diagnostic executables
# Apply SIMPLE-specific compile options to diagnostic executables
target_compile_options(diag_meiss.x PRIVATE ${SIMPLE_COMPILE_OPTIONS})
target_compile_options(diag_albert.x PRIVATE ${SIMPLE_COMPILE_OPTIONS})
target_compile_options(diag_orbit.x PRIVATE ${SIMPLE_COMPILE_OPTIONS})

# Apply trampoline error flags only to SIMPLE executable (not subprojects like fortplot)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(diag_meiss.x PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wtrampolines>
        $<$<COMPILE_LANGUAGE:Fortran>:-Werror=trampolines>
    )
    target_compile_options(diag_albert.x PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wtrampolines>
        $<$<COMPILE_LANGUAGE:Fortran>:-Werror=trampolines>
    )
    target_compile_options(diag_orbit.x PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wtrampolines>
        $<$<COMPILE_LANGUAGE:Fortran>:-Werror=trampolines>
    )
endif()

target_link_libraries(diag_meiss.x PRIVATE simple)
target_link_libraries(diag_albert.x PRIVATE simple)
target_link_libraries(diag_newton.x PRIVATE simple)
target_link_libraries(diag_orbit.x PRIVATE simple)
set_target_properties(diag_meiss.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set_target_properties(diag_albert.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set_target_properties(diag_newton.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set_target_properties(diag_orbit.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Ensure fortplot is built before diagnostics (when available)
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    add_dependencies(diag_meiss.x fortplot)
    add_dependencies(diag_newton.x fortplot)
    add_dependencies(diag_albert.x fortplot)
    add_dependencies(diag_orbit.x fortplot)
endif()

if (ENABLE_PYTHON_INTERFACE)
	add_subdirectory(python)
endif ()

if (SIMPLE_TESTING)
	enable_testing()
    add_subdirectory(test)
endif ()
