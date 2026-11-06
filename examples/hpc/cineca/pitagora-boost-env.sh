#!/bin/bash
# Source from the SIMPLE root on the Pitagora Booster, in an initialized module shell.
module purge
module load profile/base gcc/12.3.0 cuda/12.6 cmake/3.27.9 ninja/1.11.1 \
    openblas/0.3.26--gcc--12.3.0 python/3.11.7 netcdf-c/4.9.2--gcc--12.3.0

export CMAKE_BUILD_PARALLEL_LEVEL=4
export OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
export OPENBLAS_NUM_THREADS=1 CTEST_PARALLEL_LEVEL=1
