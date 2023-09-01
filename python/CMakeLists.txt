cmake_minimum_required(VERSION 3.18)

# Grab Python
find_package(Python COMPONENTS Interpreter Development NumPy)

if (Python_FOUND)
    message(STATUS "Python found: ${Python_EXECUTABLE}")
    message(STATUS "Python version: ${Python_VERSION}")
    message(STATUS "Python include dirs: ${Python_INCLUDE_DIRS}")
    message(STATUS "Python libraries: ${Python_LIBRARIES}")
    message(STATUS "Python NumPy include dirs: ${Python_NumPy_INCLUDE_DIRS}")
else()
    message("Python with NumPy not found, skipping interface build.")
    return()
endif()

# Check if f90wrap package ist installed
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m "f90wrap"
    RESULT_VARIABLE F90WRAP_RESULT
    OUTPUT_QUIET
    ERROR_QUIET
)

if (F90WRAP_RESULT EQUAL 0)
    message(STATUS "Python f90wrap found.")
else()
    message(STATUS "Python f90wrap not found, skipping interface build.")
    return()
endif()

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import numpy; print(numpy.__version__)"
    OUTPUT_VARIABLE NUMPY_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if (${NUMPY_VERSION} VERSION_GREATER_EQUAL "1.21.1")
	execute_process(
        COMMAND "${Python_EXECUTABLE}"
        -c "import numpy.f2py; print(numpy.f2py.get_include())"
        OUTPUT_VARIABLE F2PY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
	)
else()
    execute_process(
        COMMAND "${Python_EXECUTABLE}"
        -c "import numpy; print(numpy.get_include())"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(F2PY_INCLUDE_DIR "${NUMPY_INCLUDE_DIR}/../../f2py/src/")
endif()

include_directories(
    BEFORE
    ${CMAKE_BINARY_DIR}/SRC
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
    ${F2PY_INCLUDE_DIR}
)

message(STATUS "Python include dir: ${Python_INCLUDE_DIRS}")
message(STATUS "Python f2py include dir: ${F2PY_INCLUDE_DIR}")
message(STATUS "Python numpy include dir: ${Python_NumPy_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python_LIBRARIES}")
message(STATUS "Python binary output dir: ${CMAKE_CURRENT_BINARY_DIR}")

set(SOURCES_TO_WRAP
    ${PROJECT_SOURCE_DIR}/SRC/params.f90
    ${PROJECT_SOURCE_DIR}/SRC/field_can.f90
    ${PROJECT_SOURCE_DIR}/SRC/orbit_symplectic.f90
    ${PROJECT_SOURCE_DIR}/SRC/get_canonical_coordinates.f90
    ${PROJECT_SOURCE_DIR}/SRC/canonical_coordinates_mod.f90
    ${PROJECT_SOURCE_DIR}/SRC/spline_vmec_data.f90
    ${PROJECT_SOURCE_DIR}/SRC/simple.f90
    ${PROJECT_SOURCE_DIR}/SRC/bench.f90
)

set(FILES_TO_WRAP
    ${CMAKE_BINARY_DIR}/params.f90.i
    ${CMAKE_BINARY_DIR}/field_can.f90.i
    ${CMAKE_BINARY_DIR}/orbit_symplectic.f90.i
    ${CMAKE_BINARY_DIR}/get_canonical_coordinates.f90.i
    ${CMAKE_BINARY_DIR}/canonical_coordinates_mod.f90.i
    ${CMAKE_BINARY_DIR}/spline_vmec_data.f90.i
    ${CMAKE_BINARY_DIR}/simple.f90.i
    ${CMAKE_BINARY_DIR}/bench.f90.i
)

add_custom_target(PREPROCESS
    DEPENDS ${SOURCES_TO_WRAP}
)

foreach (SOURCE_TO_WRAP ${SOURCES_TO_WRAP})
get_filename_component(SOURCE_NAME ${SOURCE_TO_WRAP} NAME_WE)
add_custom_command(
    TARGET PREPROCESS
    COMMAND gfortran -E -cpp
        -o ${CMAKE_BINARY_DIR}/${SOURCE_NAME}.f90.i
        ${SOURCE_TO_WRAP}
    DEPENDS ${SOURCES_TO_WRAP}
)
endforeach()

set(F90WRAP_OUTPUT
    f90wrap_simple.f90.f90
    f90wrap_orbit_symplectic.f90.f90
    f90wrap_field_can.f90.f90
    f90wrap_get_canonical_coordinates.f90.f90
    f90wrap_params.f90.f90
    f90wrap_canonical_coordinates_mod.f90.f90
    f90wrap_bench.f90.f90
    f90wrap_toplevel.f90
)

set(F2PY_MODULE_NAME "_simple")
set(F2PY_MODULE_C "${F2PY_MODULE_NAME}module.c")
set(GENERATED_MODULE_FILE "${F2PY_MODULE_NAME}.${Python_SOABI}")
set(F2PY_WRAPPERS "${CMAKE_CURRENT_BINARY_DIR}/${F2PY_MODULE_C}")

add_custom_command(
	OUTPUT ${F90WRAP_OUTPUT} simple.py .f2py_f2cmap
	COMMAND ${Python_EXECUTABLE}
        -m f90wrap
        --f90wrap
        -m simple
        ${FILES_TO_WRAP}
	DEPENDS PREPROCESS
)

add_custom_command(
    OUTPUT ${F2PY_WRAPPERS}
    COMMAND ${Python_EXECUTABLE}
        -m f90wrap
        --f2py-f90wrap
        ${F90WRAP_OUTPUT}
        -m ${F2PY_MODULE_NAME}
        --lower
    DEPENDS PREPROCESS ${F90WRAP_OUTPUT}
)

add_custom_target(
    pysimple
    DEPENDS "${F2PY_MODULE_NAME}"
)

add_library(${F2PY_MODULE_NAME} SHARED
   ${F2PY_WRAPPERS}
   ${F2PY_INCLUDE_DIR}/fortranobject.c
   ${F90WRAP_OUTPUT}
)
add_dependencies(${F2PY_MODULE_NAME} SIMPLE_LIB)

target_link_libraries(${F2PY_MODULE_NAME} PRIVATE ${Python_LIBRARIES} SIMPLE_LIB)

set_target_properties(
    ${F2PY_MODULE_NAME}
    PROPERTIES
    PREFIX ""
    OUTPUT_NAME "${GENERATED_MODULE_FILE}"
    SUFFIX ".so"
    LINKER_LANGUAGE C
)
