# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.12)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project (SIMPLE)
enable_language(C Fortran)

add_compile_options(-g -fbacktrace)

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

include_directories(${NETCDFINCLUDE_DIR})
link_directories(${NETCDFLIB_DIR})
add_link_options(${NETCDF_FLIBS})

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_or_fetch(libneo)

if (SIMPLE_TESTING)
    message(STATUS "Unit Tests enabled!")
endif()

message(STATUS "CMake build type: " ${CMAKE_BUILD_TYPE})

include_directories ($ENV{NETCDFF_INCLUDE} ${NFINC})
link_directories ($ENV{NETCDF_LIB} $ENV{NETCDFF_LIB} ${NFLIBS})

# Set SIMPLE-specific compiler flags after GVEC to avoid affecting it
# GNU Fortran compiler flags
add_compile_options(-Wall -Wextra -Wimplicit-interface -Wno-external-argument-mismatch -fPIC -fmax-errors=5)
add_compile_options(-Wno-unused-dummy-argument)
add_compile_options(-Wtrampolines)
add_compile_options(-ffree-line-length-132)
add_compile_options(-g -fbacktrace)

# TODO: We have to allow executable stack until field_can_meiss works without it.
#       See warnings from -Wtrampolines
# if (NOT APPLE)
#     add_link_options(-Wl,-z,noexecstack)
# endif()

if (ENABLE_OPENMP)
	add_compile_options(-fopenmp)
	add_link_options(-fopenmp)
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
                    '*/libneo/*'
                    '*code/libneo/*'
                    '*/../*'
                    '*/test/*'
                    '/usr/*'
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
	add_compile_options(-C)
	add_compile_options(-O0 -ggdb
		-ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all,no-array-temps
		-finit-real=nan -fno-omit-frame-pointer -fno-inline)
elseif (CMAKE_BUILD_TYPE MATCHES Profile)
	add_compile_options(-O2 -p -shared-libgcc)
elseif (CMAKE_BUILD_TYPE MATCHES Release)
	add_compile_options(-O3)
elseif (CMAKE_BUILD_TYPE MATCHES Fast)
	add_compile_options(-O3 -ffast-math -ffp-contract=fast -funroll-loops -mtune=native)
	if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64")
			add_compile_options(-march=native)
	else()
	add_compile_options(-march=haswell)
	endif()
endif()

add_subdirectory(src)
include_directories("${CMAKE_BINARY_DIR}/src")

add_executable (simple.x
	app/simple.f90
)

target_link_libraries(simple.x PRIVATE simple)
set_target_properties(simple.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (ENABLE_PYTHON_INTERFACE)
	add_subdirectory(python)
endif ()

if (SIMPLE_TESTING)
	enable_testing()
    add_subdirectory(test)
endif ()
