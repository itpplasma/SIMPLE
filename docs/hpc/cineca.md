# CINECA Leonardo and Pitagora

This page builds, tests and runs SIMPLE on the CINECA machines Leonardo and
Pitagora. Tested on Leonardo (CPU, GCC 12.2, OpenMP) at SIMPLE `fe2de27`:
build, the test selection and the smoke run below pass, also as a batch job.
The CUDA build compiles but was not run on a GPU. The Pitagora recipe has not
been run.

## Access

Follow CINECA's [access and 2FA instructions](https://docs.hpc.cineca.it/general/access.html):
enrol an OTP authenticator, install the Smallstep `step` client and bootstrap
its certificate authority with the published URL and fingerprint. Then, with
an SSH agent running:

```bash
step ssh login '<registered-email>' --provisioner cineca-hpc
ssh <username>@login.leonardo.cineca.it
```

Certificates last 12 hours. Cloning SIMPLE from GitHub over HTTPS works on
the login nodes; CMake fetches libneo, Fortnum and Fortplot the same way.

## Leonardo: build

`examples/hpc/cineca/leonardo-gcc-env.sh`, sourced from the SIMPLE root,
loads the module stack and limits threads for the shared login node:

```bash
module purge
module load profile/base gcc/12.2.0 cmake/3.27.9 ninja/1.11.1 \
    openblas/0.3.26--gcc--12.2.0 python/3.11.7 \
    netcdf-c/4.9.2--gcc--12.2.0-spack0.22
export CMAKE_BUILD_PARALLEL_LEVEL=4
export OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
export OPENBLAS_NUM_THREADS=1 CTEST_PARALLEL_LEVEL=1
export TMPDIR="$PWD/build/tmp"; mkdir -p "$TMPDIR"
```

MPI, NetCDF-Fortran, FFTW, GSL and MKL are not needed: NetCDF-C (which loads
serial HDF5) serves SIMPLE's C test-data generator, and Fortran NetCDF access
goes through Fortio, fetched by libneo. Python 3 is needed for the test
drivers even with the Python interface off.

```bash
git clone https://github.com/itpplasma/SIMPLE.git && cd SIMPLE
source examples/hpc/cineca/leonardo-gcc-env.sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_OPENMP=ON -DSIMPLE_DETERMINISTIC_FP=ON \
    -DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF \
    -DLIBNEO_REF=f6090e8150c2f5b1c0afde09baadf79b50b1ffe7
cmake --build build -j 4
```

Use Release, not `Fast`, which tunes for the build host. Without
`LIBNEO_REF` CMake takes the libneo branch matching SIMPLE's branch name, or
`main`; the value above is the tested pin. The clean build took 107 s.

## Leonardo: test and run

On the login node, run a short selection of the fast suite (11 tests, 13 s),
keeping the libneo pin in `FLAGS`:

```bash
make test-fast CONFIG=Release \
    FLAGS='-DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF -DENABLE_OPENMP=ON -DLIBNEO_REF=f6090e8150c2f5b1c0afde09baadf79b50b1ffe7' \
    TEST='^(test_config_reference|test_finite_check|test_profiles|test_deterministic_seed|test_orbit_invariants|test_stencil_utils|test_orbit_symplectic_base|test_newton_solver_status|test_spectre_fo_spline_tok2vol|test_spectre_fo_spline_ncsx_freeboundary|test_diag_albert)$'
```

The full fast suite contains integration cases with 900 s timeouts; run it in
a batch job, not on the login node.

The smoke run traces eight markers for $10^{-5}$ s in the QA equilibrium that
the test configuration downloads to `test/test_data/wout.nc`:

```bash
mkdir -p run/cineca-smoke
cp examples/hpc/cineca/smoke.in run/cineca-smoke/simple.in
cp test/test_data/wout.nc run/cineca-smoke/wout.nc
(cd run/cineca-smoke && ../../build/simple.x simple.in)
```

Expected: exit status 0 after about 10 s, 8 rows in `times_lost.dat` and
`orbit_exit_code.dat`, 20 rows in `confined_fraction.dat`, all finite, final
time `1e-5`; exit codes four `0` (traced) and four `3` (skipped as strongly
passing, default `contr_pp`).

As a batch job on the budget-free serial partition (no account needed; one
task, four cores, 4 GB, three minutes), submitted from the SIMPLE root:

```bash
sbatch examples/hpc/cineca/leonardo-serial.sbatch
sacct -j <job-id> --format=JobID,State,ExitCode,Elapsed
```

Production jobs follow {doc}`index`; partitions and limits are in
[Leonardo's Slurm guide](https://docs.hpc.cineca.it/hpc/leonardo.html#job-managing-and-slurm-partitions).

## Leonardo: CUDA

Run GPU jobs only in an allocation with a GPU, never on login or serial
nodes. Configuration and valid inputs: {doc}`/user-guide/gpu`. The A100 build
uses the CPU stack plus CUDA 12.3:

```bash
module load cuda/12.3
cmake -S . -B build-cuda -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DSIMPLE_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF \
    -DLIBNEO_REF=f6090e8150c2f5b1c0afde09baadf79b50b1ffe7
cmake --build build-cuda -j 4 --target \
    simple.x test_cuda_native_rk54.x test_cuda_native_bridge.x
cuobjdump --list-elf build-cuda/libsimple_cuda_native_rk54.so   # expect sm_80
```

The device code lives in the shared library `libsimple_cuda_native_rk54.so`.
Validate with `make test ... TEST='test_cuda_native_(rk54|bridge)'` inside a
GPU job before production.

NVHPC 24.5 (`nvfortran`) does not build SIMPLE: Fortnum's
`src/verified/fortnum_interval_qp.f90` declares `real(real128)`, for which
NVFORTRAN returns kind `-3`. Use GNU Fortran with `nvcc`.

## Leonardo: pitfalls

- `netcdf.h` not found in `test/tests/generate_test_boozer_chartmap.c`: load
  the `netcdf-c` module.
- `module load git/2.38.1` together with the OpenMPI stack fails (its
  `openssh` dependency loads an older Kerberos). Use `/usr/bin/git` after
  `module purge`.

## Pitagora (untested)

Consult the current [Pitagora guide](https://docs.hpc.cineca.it/hpc/pitagora.html).
Use `login-dcgp.pitagora.cineca.it` for the AMD Turin CPU partition and
`login-boost.pitagora.cineca.it` for the Intel/H100 Booster. Leonardo's
module names are not known to exist there: find GCC, CMake, Ninja, Python,
OpenBLAS and NetCDF-C with `modmap -m <name>`, write the loads into an
environment script like `leonardo-gcc-env.sh`, export its absolute path as
`SIMPLE_ENV`, and repeat the Leonardo build, test selection and smoke run on
the destination partition. For the H100 use `-DCMAKE_CUDA_ARCHITECTURES=90`.

| Partition | Nodes | Max. nodes/job | Max. wall time |
| --- | --- | ---: | --- |
| `dcgp_fua_dbg` / `dcgp_fua_prod` | 2 × AMD EPYC 9745, 256 cores, 768 GiB | 2 / 64 | 30 min / 24 h |
| `boost_fua_dbg` / `boost_fua_prod` | 64 cores, 4 × H100, 512 GiB | 2 / 16 | 30 min / 24 h |
| `ptgr_all_serial` | CPU only, budget-free, max. 4 cores | — | 4 h |

One SIMPLE process uses one node or one GPU. Templates (replace the account
before submission):

| Template | Request | Submit from |
| --- | --- | --- |
| `examples/hpc/cineca/pitagora-cpu.sbatch` | debug, 4 cores, 3 min | SIMPLE root, runs `run/cineca-smoke` |
| `examples/hpc/cineca/pitagora-gpu.sbatch` | debug, 1 GPU, 3 min | SIMPLE root, runs `run/cineca-gpu` |
| `examples/hpc/pitagora/dcgp.slurm` | debug, 256 cores, 30 min | run directory; `SIMPLE_EXE` = path to `simple.x` |
| `examples/hpc/pitagora/booster.slurm` | debug, 1 H100 + 16 cores, 30 min | run directory; `SIMPLE_EXE` = path to `build-cuda/simple.x` |

For production, override on the command line, e.g.
`sbatch --partition=dcgp_fua_prod --time=24:00:00 examples/hpc/pitagora/dcgp.slurm`.
A GPU input sets `gpu_mode = 'production'`, `gpu_backend = 'cuda-native'`,
`gpu_method = 'dopri'`, `gpu_num_devices = 1`, `isw_field_type = 2`, with
collisions, walls, classifiers and orbit output off. Native CUDA tracing
needs a Boozer chartmap field; with a VMEC file it stops before tracing
(`native CUDA RK requires initialized compact Boozer tables`). Compare seeded
CPU and GPU results on the same case before production.
