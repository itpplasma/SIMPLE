
set (SOURCES
    # RKF45 from https://people.sc.fsu.edu/~jburkardt/
    # licenced under LGPL (see LICENSE.rkf45)
    minpack.f90
)

add_library(RKF45 SHARED
    rkf45.f90
)

if ($ENV{SIMPLE_COVERAGE} STREQUAL "TRUE")
    add_subdirectory(pFUnit)
endif()

add_library(CONTRIB SHARED ${SOURCES})
target_link_libraries(CONTRIB RKF45)

if ($ENV{SIMPLE_COVERAGE} STREQUAL "TRUE")
    target_link_libraries(CONTRIB pfunit)
endif()



