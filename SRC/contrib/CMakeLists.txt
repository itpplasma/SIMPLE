
set (SOURCES
    # RKF45 from https://people.sc.fsu.edu/~jburkardt/
    # licenced under LGPL (see LICENSE.rkf45)
    minpack.f90
)

add_library(RKF45 SHARED
    rkf45.f90
)

add_library(CONTRIB SHARED ${SOURCES})
target_link_libraries(CONTRIB RKF45)




