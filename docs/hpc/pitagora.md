# Pitagora deployment templates (untested)

These recipes and [CPU](../../examples/hpc/pitagora/dcgp.slurm) and
[GPU](../../examples/hpc/pitagora/booster.slurm) jobs have **not been tested on
Pitagora**. Local workstation results are in
[the scaling report](../performance/efss-hpc/index.rst).

## Public allocation information

CINECA's [Pitagora guide](https://docs.hpc.cineca.it/hpc/pitagora.html),
consulted 2026-10-18, documents the following normal-QOS limits:

| Partition | Maximum nodes/job | Maximum wall time |
| --- | ---: | ---: |
| `dcgp_fua_dbg` | 2 | 00:30:00 |
| `dcgp_fua_prod` | 64 | 24:00:00 |
| `boost_fua_dbg` | 2 | 00:30:00 |
| `boost_fua_prod` | 16 | 24:00:00 |

A DCGP node has two AMD Turin EPYC 9745 processors, 256 cores and 768 GiB RAM.
A Booster node has 64 CPU cores, four NVIDIA H100 GPUs and 512 GiB host RAM.
These are published specifications, not SIMPLE measurements. The templates
request one node, one process and either 256 CPU cores or one GPU plus 16 cores.
OpenMP shares memory within a node; these jobs do not distribute one simulation
across nodes. Use separate run directories for independent ensemble jobs.

## Toolchain and build

Use `login-dcgp.pitagora.cineca.it` for CPU builds and
`login-boost.pitagora.cineca.it` for GPU builds. Check module availability first:

```bash
module purge
modmap -m gcc
modmap -m netcdf-fortran
modmap -m openblas
modmap -m cmake
modmap -m ninja
modmap -m cuda
```

CINECA's [module guide](https://docs.hpc.cineca.it/hpc/hpc_environment.html)
describes `modmap` and compiler selection. A published
[module inventory](https://docs.hpc.cineca.it/_downloads/84493ea9f614a4a9b4d303f0fb5f34cd/Pitagora_EF.pdf)
lists GCC 12.3.0, NetCDF-C 4.9.2, NetCDF-Fortran 4.6.1 and CMake 3.27.9.
That inventory is historical; **resolve compatible installed versions before
using the following template**. Do not mix GNU and NVHPC Fortran module files.

```bash
module load gcc/<version> cmake/<version> ninja/<version>
module load netcdf-c/<gcc-compatible-version> netcdf-fortran/<gcc-compatible-version>
module load openblas/<gcc-compatible-version>
export CC=gcc CXX=g++ FC=gfortran
export CMAKE_BUILD_PARALLEL_LEVEL=16 OPENBLAS_NUM_THREADS=1
# All commands below run from the SIMPLE repository root.
make CONFIG=Release FLAGS='-DENABLE_OPENMP=ON -DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF'
make test-fast CONFIG=Release FLAGS='-DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF'
module -t list 2> build/modules.txt
gfortran --version > build/compiler.txt
git rev-parse HEAD > build/source-commit.txt
```

Dependencies fetched by CMake require network access during configuration;
pre-stage them where permitted. Preserve their actual commits and `CMakeCache.txt`
alongside the executable. SIMPLE needs NetCDF, BLAS/LAPACK and dependencies
resolved by CMake; let the configuration fail visibly if the module stack is incomplete.
The test target enables deterministic floating point; rebuild with the intended
`SIMPLE_DETERMINISTIC_FP` setting before comparing performance across builds.

For the native CUDA path, keep the same GNU stack and add a compatible CUDA module:

```bash
module load cuda/<gcc-compatible-version>
cmake -S . -B build-cuda -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DSIMPLE_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 \
  -DENABLE_OPENMP=ON -DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF
cmake --build build-cuda --parallel 16
```

H100 uses architecture `90`. Native CUDA uses a C++17 CUDA compiler plus GNU
Fortran; NVIDIA Fortran is optional. Validate both native CUDA tests before production use
(the Makefile's command-line `BUILD_DIR` override keeps its native test driver):

```bash
OMP_NUM_THREADS=16 CMAKE_BUILD_PARALLEL_LEVEL=16 \
  make test CONFIG=Release BUILD_DIR=build-cuda \
  FLAGS='-DSIMPLE_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 -DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF' \
  TEST='test_cuda_native_(rk54|bridge)'
```

Run this inside an allocated GPU job. It enables deterministic floating point;
reconfigure the intended setting before timing. See
[GPU controls](../../DOC/gpu-cuda.md). Use Boozer guiding-center tracing,
`gpu_mode='production'`, `gpu_backend='cuda-native'`, `gpu_method='dopri'`,
`gpu_landreman=.true.`, and `gpu_num_devices=1`. The local direct-VMEC production attempt currently
fails at `src/simple_cuda_native.f90:133` before CUDA tracing; see the linked
performance report. This template does not resolve that initialization block.
Validate a compatible chartmap/Boozer-start workflow or an upstream fix before
production deployment. Collisions, walls, classifiers
and orbit output are unsupported. The ordinary CPU template uses `gpu_mode='off'`.

## Submit and validate

Replace `--account=<project>` in your copy with the granted project account.
Prepare a writable run directory with `simple.in` and its equilibrium, load the
build's module stack, and submit there:

```bash
export SIMPLE_EXE=/absolute/path/to/SIMPLE/build/simple.x
sbatch /absolute/path/to/SIMPLE/examples/hpc/pitagora/dcgp.slurm
# Production alternative; command-line options override the debug headers:
sbatch --partition=dcgp_fua_prod --time=24:00:00 /absolute/path/to/SIMPLE/examples/hpc/pitagora/dcgp.slurm
# GPU alternative:
export SIMPLE_EXE=/absolute/path/to/SIMPLE/build-cuda/simple.x
sbatch /absolute/path/to/SIMPLE/examples/hpc/pitagora/booster.slurm
sbatch --partition=boost_fua_prod --time=24:00:00 /absolute/path/to/SIMPLE/examples/hpc/pitagora/booster.slurm
```

Start in debug. Check the allocation, thread binding, selected CUDA backend,
`confined_fraction.dat`, `unresolved_fraction.dat`, and `orbit_exit_code.dat`.
Compare seeded CPU/GPU results on the same public case before a production run.
Then switch to the corresponding production partition and an appropriate wall
time within the limits above. `bash -n` checks local shell syntax only; it does
not validate Slurm resources, module availability or numerical execution.

Open questions requiring access: actual module versions, NetCDF/BLAS discovery,
Slurm GPU binding, memory consumption and scaling across the two CPU sockets.
No Pitagora throughput or efficiency claim is made by these templates.
