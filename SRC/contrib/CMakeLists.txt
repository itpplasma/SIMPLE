
set (SOURCES
    # RKF45 from https://people.sc.fsu.edu/~jburkardt/
    # licenced under LGPL (see LICENSE.rkf45)
    rkf45.f90
    minpack.f90
)

if ($ENV{SIMPLE_COVERAGE} STREQUAL "TRUE")
    add_subdirectory(pFUnit)
endif()

add_library(CONTRIB PRIVATE ${SOURCES})

if ($ENV{SIMPLE_COVERAGE} STREQUAL "TRUE")
    target_link_libraries(CONTRIB pFUnit)
endif()

