
set (SOURCES
    # RKF45 from https://people.sc.fsu.edu/~jburkardt/
    # licenced under LGPL (see LICENSE.rkf45)
    minpack_interfaces.f90
    minpack.f90
)
add_library(CONTRIB STATIC ${SOURCES})

# Disable specific warnings for minpack.f90 as it's legacy code
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_source_files_properties(minpack.f90 PROPERTIES
        COMPILE_FLAGS "-Wno-compare-reals -Wno-implicit-interface")
endif()
