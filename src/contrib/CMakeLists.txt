
set (SOURCES
    # RKF45 from https://people.sc.fsu.edu/~jburkardt/
    # licenced under LGPL (see LICENSE.rkf45)
    minpack.f90
)
add_library(CONTRIB STATIC ${SOURCES})



