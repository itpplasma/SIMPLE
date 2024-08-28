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
    cut_detector.f90
    orbit_symplectic.f90
    field_can.f90
    get_canonical_coordinates.f90
    params.f90
    canonical_coordinates_mod.f90
    vmecinm_m.f90
    new_vmec_allocation_stuff.f90
    spline_vmec_data.f90
    magfie.f90
)

set(EXCLUDE_BASENAMES
    # Empty
)

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

list(TRANSFORM FILES_TO_WRAP PREPEND ${PROJECT_SOURCE_DIR}/SRC/)

set(F2PY_MODULE_NAME "_pysimple")

set(F2PY_MODULE_C "${F2PY_MODULE_NAME}module.c")
set(GENERATED_MODULE_FILE "${F2PY_MODULE_NAME}.${Python_SOABI}")
set(F2PY_WRAPPERS "${CMAKE_CURRENT_BINARY_DIR}/${F2PY_MODULE_C}")

add_custom_target(F90WrapOutput)

foreach (file_to_wrap ${FILES_TO_WRAP})
    get_filename_component(basename ${file_to_wrap} NAME_WE)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${basename}.f90.i
        COMMAND gfortran -E -cpp
            -o ${CMAKE_CURRENT_BINARY_DIR}/${basename}.f90.i
            ${file_to_wrap}
    DEPENDS ${file_to_wrap}
)
endforeach()

add_custom_command(OUTPUT ${F90WRAP_OUTPUT}
	COMMAND ${Python_EXECUTABLE}
        -m f90wrap
        --f90wrap
        -m pysimple
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
    OUTPUT ${CMAKE_BINARY_DIR}/pysimple.py
    DEPENDS "${F2PY_MODULE_NAME}"
    COMMENT "copying file"
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/pysimple.py
            ${CMAKE_BINARY_DIR}/pysimple.py
)

add_custom_target(
    pysimple ALL
    DEPENDS ${CMAKE_BINARY_DIR}/pysimple.py
)
