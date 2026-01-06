cmake_minimum_required(VERSION 3.22)

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
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
    ${F2PY_INCLUDE_DIR}
)

message(STATUS "Python include dir: ${Python_INCLUDE_DIRS}")
message(STATUS "Python f2py include dir: ${F2PY_INCLUDE_DIR}")
message(STATUS "Python numpy include dir: ${Python_NumPy_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python_LIBRARIES}")
message(STATUS "Python binary output dir: ${CMAKE_CURRENT_BINARY_DIR}")

set(FILES_TO_WRAP
    simple_main.f90
    simple.f90
    samplers.f90
    cut_detector.f90
    callback.f90
    orbit_symplectic_base.f90
    orbit_symplectic.f90
    field_can.f90
    get_canonical_coordinates.F90
    coordinates/stencil_utils.f90
    coordinates/array_utils.f90
    params.f90
    params_wrapper.f90
    magfie_wrapper.f90
)

# Add libneo files to wrap
set(LIBNEO_SOURCE_DIR ${libneo_SOURCE_DIR})

list(APPEND FILES_TO_WRAP ${LIBNEO_SOURCE_DIR}/src/canonical_coordinates_mod.f90)

set(EXCLUDE_BASENAMES
    # Empty
)

# Separate SIMPLE and libneo files
set(SIMPLE_FILES_TO_WRAP)
set(LIBNEO_FILES_TO_WRAP)

foreach(file ${FILES_TO_WRAP})
    if(file MATCHES "^${LIBNEO_SOURCE_DIR}")
        list(APPEND LIBNEO_FILES_TO_WRAP ${file})
    else()
        list(APPEND SIMPLE_FILES_TO_WRAP ${file})
    endif()
endforeach()

# Add PROJECT_SOURCE_DIR/src/ prefix to SIMPLE files only
list(TRANSFORM SIMPLE_FILES_TO_WRAP PREPEND ${PROJECT_SOURCE_DIR}/src/)

# Combine the lists
set(FILES_TO_WRAP ${SIMPLE_FILES_TO_WRAP} ${LIBNEO_FILES_TO_WRAP})

foreach(file ${FILES_TO_WRAP})
    get_filename_component(basename ${file} NAME_WE)
    list(APPEND PREPROCESSED_SOURCES
        ${CMAKE_CURRENT_BINARY_DIR}/${basename}.f90.i
    )
    if(NOT ${basename} IN_LIST EXCLUDE_BASENAMES)
        list(APPEND F90WRAP_OUTPUT
            ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_${basename}.f90.f90
        )
    endif()
endforeach()

set(F2PY_MODULE_NAME "_simple_backend")

set(F2PY_MODULE_C "${F2PY_MODULE_NAME}module.c")
set(GENERATED_MODULE_FILE "${F2PY_MODULE_NAME}.${Python_SOABI}")
set(F2PY_WRAPPERS "${CMAKE_CURRENT_BINARY_DIR}/${F2PY_MODULE_C}")

set(SIMPLE_F90PP_DEFINES)
if(ENABLE_GVEC)
    list(APPEND SIMPLE_F90PP_DEFINES -DGVEC_AVAILABLE)
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    list(APPEND SIMPLE_F90PP_DEFINES -DSIMPLE_NVHPC)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(SIMPLE_F90PP_FLAGS -Mpreprocess -E)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(SIMPLE_F90PP_FLAGS -E -cpp)
else()
    set(SIMPLE_F90PP_FLAGS -E)
endif()

string(JOIN " " SIMPLE_F90PP_FLAGS_STR ${SIMPLE_F90PP_FLAGS})
string(JOIN " " SIMPLE_F90PP_DEFINES_STR ${SIMPLE_F90PP_DEFINES})

foreach (file_to_wrap ${FILES_TO_WRAP})
    get_filename_component(basename ${file_to_wrap} NAME_WE)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${basename}.f90.i
        COMMAND bash -c
            "${CMAKE_Fortran_COMPILER} ${SIMPLE_F90PP_FLAGS_STR} ${SIMPLE_F90PP_DEFINES_STR} \"${file_to_wrap}\" > \"${CMAKE_CURRENT_BINARY_DIR}/${basename}.f90.i\""
        DEPENDS ${file_to_wrap}
        VERBATIM
    )
endforeach()

add_custom_command(OUTPUT ${F90WRAP_OUTPUT}
	COMMAND ${Python_EXECUTABLE}
        -m f90wrap
        --f90wrap
        -m simple_backend
        -k ${CMAKE_CURRENT_SOURCE_DIR}/f2py_f2cmap
        ${CMAKE_CURRENT_BINARY_DIR}/*.f90.i
        > ${CMAKE_CURRENT_BINARY_DIR}/f90wrap.log 2>&1
	DEPENDS ${PREPROCESSED_SOURCES}
    COMMENT "Processing preprocessed sources with f90wrap"
)

add_custom_command(OUTPUT ${F2PY_WRAPPERS}
    COMMAND ${Python_EXECUTABLE}
        -m f90wrap
        --f2py-f90wrap
        ${F90WRAP_OUTPUT}
        -m ${F2PY_MODULE_NAME}
        --f2cmap ${CMAKE_CURRENT_SOURCE_DIR}/f2py_f2cmap
        --lower
        > ${CMAKE_CURRENT_BINARY_DIR}/f2py.log 2>&1
    DEPENDS ${F90WRAP_OUTPUT}
    COMMENT "Processing wrapped sources with f2py"
)

add_library(${F2PY_MODULE_NAME} SHARED
   ${F2PY_WRAPPERS}
   ${F2PY_INCLUDE_DIR}/fortranobject.c
   ${F90WRAP_OUTPUT}
)

add_dependencies(${F2PY_MODULE_NAME} simple)

target_link_libraries(${F2PY_MODULE_NAME} PRIVATE ${Python_LIBRARIES} simple)

set_target_properties(
    ${F2PY_MODULE_NAME}
    PROPERTIES
    PREFIX ""
    OUTPUT_NAME "${GENERATED_MODULE_FILE}"
    SUFFIX ".so"
    LINKER_LANGUAGE C
)

add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/simple_backend.py
    DEPENDS "${F2PY_MODULE_NAME}"
    COMMENT "copying file"
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/simple_backend.py
            ${CMAKE_BINARY_DIR}/simple_backend.py
)

add_custom_target(
    simple_backend ALL
    DEPENDS ${CMAKE_BINARY_DIR}/simple_backend.py
)

install(TARGETS ${F2PY_MODULE_NAME} DESTINATION .)
install(FILES ${CMAKE_BINARY_DIR}/simple_backend.py DESTINATION .)
