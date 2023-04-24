# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.12)
project (SIMPLE Fortran C)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
#    set(CMAKE_BUILD_TYPE Debug)
endif()

option(ENABLE_PYTHON_INTERFACE "Enables the Python-Wrapper." OFF)
option(SIMPLE_TESTING "Enable testing." OFF)

set(CMAKE_MACOSX_RPATH 1)

find_program(NF_CONFIG "nc-config")

if (NF_CONFIG)
    execute_process(COMMAND "nc-config" --prefix
		OUTPUT_VARIABLE NFPREFIX)
else()
    message(SEND_ERROR "nc-config not found. Please install libnetcdff-dev")
endif()

string(STRIP ${NFPREFIX} NFPREFIX)
set(NFINC ${NFPREFIX}/include)
set(NFLIBS ${NFPREFIX}/lib)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)

if (SIMPLE_TESTING)
    message(STATUS "Unit Tests and Code Coverage Enabled!")
	find_package(PFUNIT REQUIRED)
	enable_testing()
	add_link_options(--coverage)
endif()

message(STATUS "CMake build type: " ${CMAKE_BUILD_TYPE})
message(STATUS "NetCDF include path: " ${NFINC})
message(STATUS "NetCDF lib path: " ${NFLIBS})
message(STATUS "MPI include path: " ${MPI_Fortran_INCLUDE_DIRS})

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
include_directories ($ENV{NETCDFF_INCLUDE} ${NFINC} ${MPI_Fortran_INCLUDE_DIRS})
link_directories ($ENV{NETCDF_LIB} $ENV{NETCDFF_LIB} ${NFLIBS} $ENV{HOME}/.local/lib)

add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-cpp>)
add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fopenmp>)
#add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008>)
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	add_compile_options(-mkl -qopenmp -warn all,nounused -stand f2008)
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
	add_compile_options(-Wall -Wno-unused -Wno-unused-dummy-argument)
	add_link_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fopenmp>)
	add_link_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-cpp>)
	if (CMAKE_BUILD_TYPE MATCHES Debug)
		add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-C>)
		add_compile_options(-O0 -g -ggdb -fbacktrace
			-ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all,no-array-temps)
	elseif (CMAKE_BUILD_TYPE MATCHES Profile)
		add_compile_options(-O2 -p -g -shared-libgcc)
	elseif (CMAKE_BUILD_TYPE MATCHES Release)
		add_compile_options(-O3 -g -fbacktrace)
	elseif (CMAKE_BUILD_TYPE MATCHES Fast)
		add_compile_options(-O3 -march=native -mtune=native)
	endif()
endif()

add_subdirectory(SRC)

add_library (SIMPLE_LIB INTERFACE)
target_compile_options(SIMPLE_LIB INTERFACE -DMPI)
target_link_libraries(SIMPLE_LIB INTERFACE SIMPLE_SRC)
target_include_directories(SIMPLE_LIB INTERFACE ${CMAKE_BINARY_DIR}/mod/) #Why?

 #add_library (simple_mpi SHARED)
 #target_compile_options(simple_mpi PUBLIC -DMPI)
 #target_sources (simple_mpi PRIVATE ${simple_sources})
 #set_property(TARGET simple_mpi PROPERTY Fortran_MODULE_DIRECTORY mod_mpi/)
 #target_include_directories(simple_mpi PRIVATE ${CMAKE_BINARY_DIR}/mod_mpi/)
 #target_link_libraries(simple_mpi ${MPI_Fortran_LIBRARIES})

 #add_executable (simple_mpi.x
 #	main.f90
 #)
 #target_link_libraries(simple_mpi.x simple_mpi)
 #get_target_property(SIMPLE_MPI_LINKED_LIBS simple_mpi.x LINK_LIBRARIES)
 #list(REMOVE_ITEM SIMPLE_MPI_LINKED_LIBS SIMPLE_LIB)
 #set_property(TARGET simple_mpi.x PROPERTY LINK_LIBRARIES ${SIMPLE_MPI_LINKED_LIBS})

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	target_link_libraries(SIMPLE_LIB PRIVATE netcdf netcdff)
else()
	target_link_libraries(SIMPLE_LIB INTERFACE
		netcdf netcdff ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
	)
endif()

add_executable (simple.x
	SRC/main.f90
)

target_link_libraries(simple.x PRIVATE SIMPLE_LIB)

if (SIMPLE_TESTING)
    add_subdirectory(test)
endif ()


