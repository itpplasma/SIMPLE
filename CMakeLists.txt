# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.12)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project (SIMPLE)
enable_language(C Fortran)

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
link_directories ($ENV{NETCDF_LIB} $ENV{NETCDFF_LIB} ${NFLIBS} $ENV{HOME}/.local/lib)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	add_compile_options(-mkl -qopenmp -warn all,nounused)
	link_libraries("-mkl -qopenmp")
	if (CMAKE_BUILD_TYPE MATCHES Debug)
		add_compile_options(-O0 -g -traceback)
	elseif (CMAKE_BUILD_TYPE MATCHES Profile)
		add_compile_options(-O2 -g -shared-intel -debug inline-debug-info
			-D TBB_USE_THREADING_TOOLS -qopenmp-link dynamic -parallel-source-info=2)
	elseif (CMAKE_BUILD_TYPE MATCHES Release)
		add_compile_options(-O3 -g -traceback)
	elseif (CMAKE_BUILD_TYPE MATCHES Fast)
		add_compile_options(-O3 -march=native -mtune=native)
	endif()
else()
	add_compile_options(-Wall -Wextra -Wimplicit-interface -Wno-external-argument-mismatch -fPIC -fmax-errors=5)
	add_compile_options(-Wno-unused-dummy-argument)
	# If not apple
	if (NOT APPLE)
		add_link_options(-Wl,-z,noexecstack)
	endif()
	if (ENABLE_OPENMP)
		add_compile_options(-fopenmp)
		add_link_options(-fopenmp)
	endif()
	if (CMAKE_BUILD_TYPE MATCHES Debug)
		add_compile_options(-C)
		add_compile_options(-O0 -g -ggdb
			-ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all,no-array-temps
			-finit-real=nan -fno-omit-frame-pointer -fno-inline)
			add_compile_options(-fbacktrace)
	elseif (CMAKE_BUILD_TYPE MATCHES Profile)
		add_compile_options(-O2 -p -g -shared-libgcc)
	elseif (CMAKE_BUILD_TYPE MATCHES Release)
		add_compile_options(-O3 -g)
		add_compile_options(-fbacktrace)
	elseif (CMAKE_BUILD_TYPE MATCHES Fast)
		add_compile_options(-O3 -ffast-math -ffp-contract=fast -funroll-loops -mtune=native)
		if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64")
				add_compile_options(-march=native)
		else()
		add_compile_options(-march=haswell)
		endif()
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
