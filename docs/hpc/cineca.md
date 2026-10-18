# SIMPLE at CINECA: Leonardo and Pitagora

This deployment record supports EUROfusion Standard Software criterion 1.4.
Leonardo build and test commands were checked on 2026-10-18 with SIMPLE
`fe2de273681642ed1f651d3c14c5e23e398fb458`. Pitagora instructions are derived
from that test and CINECA documentation; they have **not been run on Pitagora**.
This is installation and smoke-test evidence, not a scaling or physics validation.

## Access and resource limits

Follow CINECA's [access and 2FA instructions](https://docs.hpc.cineca.it/general/access.html)
to enroll an OTP authenticator and install the Smallstep `step` client locally.
Bootstrap its certificate authority using the URL and fingerprint published
there (`step ca bootstrap`); check the current fingerprint at the source.
Use an existing SSH agent, or start one with `eval "$(ssh-agent -s)"`, then:

```bash
step ssh login '<registered-email>' --provisioner cineca-hpc
ssh <username>@login.leonardo.cineca.it
```

Complete the browser login and OTP challenge. Certificates last 12 hours;
renew with `step ssh login` when needed. Keep credentials, private keys and
certificates outside the repository. HTTPS cloning from GitHub worked from
Leonardo; no GitHub SSH key or forwarded agent is needed for the build.

These tests used four build jobs, at most four OpenMP threads, one BLAS thread,
and sequential CTest execution on a shared login node. Keep login tests short;
use Slurm for production. Store builds under `$HOME` (50 GB quota in this test)
or `$SCRATCH`; remove surplus build trees and preserve source revisions/logs.
The sole job submitted here used the budget-free `lrd_all_serial` partition
without an account. No GPU or production-partition job was submitted.

## Verified Leonardo CPU recipe

From a clean module shell, the required explicit modules were:

```bash
module purge
module load profile/base gcc/12.2.0 cmake/3.27.9 ninja/1.11.1 \
    openblas/0.3.26--gcc--12.2.0 python/3.11.7 \
    netcdf-c/4.9.2--gcc--12.2.0-spack0.22
export CMAKE_BUILD_PARALLEL_LEVEL=4
export OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
export OPENBLAS_NUM_THREADS=1 CTEST_PARALLEL_LEVEL=1
```

The reusable equivalent is
[`leonardo-gcc-env.sh`](../../examples/hpc/cineca/leonardo-gcc-env.sh), sourced
from the SIMPLE root. NetCDF-C automatically loads serial HDF5 1.14.3 and its
compression dependencies. MPI, NetCDF-Fortran, FFTW, GSL, MKL and Anaconda are
not required for this pinned native build. Their availability was inspected;
an initial build with the larger GCC/OpenMPI/NetCDF-Fortran stack also passed.

```bash
git clone https://github.com/itpplasma/SIMPLE.git "$HOME/SIMPLE-efss-cineca"
cd "$HOME/SIMPLE-efss-cineca"
git checkout --detach fe2de273681642ed1f651d3c14c5e23e398fb458
export TMPDIR="$PWD/build/tmp"
mkdir -p "$TMPDIR"
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_OPENMP=ON -DSIMPLE_DETERMINISTIC_FP=ON \
    -DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF \
    -DLIBNEO_REF=f6090e8150c2f5b1c0afde09baadf79b50b1ffe7
cmake --build build -j 4
```

The environment and input files in this guide are supplied by the documentation
branch/PR; copy them into a checkout of the tested base, or use the merged
documentation revision. Build configuration and physics sources were unchanged.
Release deliberately avoids the host-specific `Fast` profile. Deterministic
floating-point options match the repository's `make test-fast` prerequisite.

CMake fetches libneo, Fortnum and Fortplot over HTTPS. Pin libneo explicitly:
`cmake/Util.cmake:80–103` otherwise chooses a matching branch or `main`.
The resolved revisions were:

| Dependency | Revision | Selection |
| --- | --- | --- |
| libneo | `f6090e8150c2f5b1c0afde09baadf79b50b1ffe7` | Explicit recipe pin |
| Fortnum | `1a52e13c5e28b09d083947a6ec7e8a1dcf9d502e` | `cmake/dependency_versions.json` |
| Fortplot | `5daec54649e305fd09bc1ef8345511d5692a1889` | Same manifest |
| Fortio | `38add03dc691a7cd4715cd3e611ad42422e46091` | libneo's `fpm.toml`, read by its CMake |

This libneo revision fetches Fortio at `CMakeLists.txt:231–239`; SIMPLE consumes
the `fortio::fortio` target at `CMakeLists.txt:291–292`. NetCDF-Fortran and
compiler-specific HDF5 Fortran modules are consequently unnecessary here.
NetCDF-C is still required by SIMPLE's C test-data generator.

Python 3 is required for test drivers even with both Python options disabled
(`test/tests/CMakeLists.txt:47`). This recipe does not build or validate
`pysimple`; its separate prerequisites are in `requirements.txt` and
[the Python guide](../../python/README.rst). No Python packages were installed
on Leonardo for this native smoke test.

## Short tests and example

Run tests through the repository's Make targets, from its root. This bounded
selection of the fast suite covers configuration, numerical checks, native
diagnostics, and four-thread SPECTRE spline tests with real HDF5 fixtures:

```bash
/usr/bin/time -p make test-fast CONFIG=Release \
    LIBNEO_REF=f6090e8150c2f5b1c0afde09baadf79b50b1ffe7 \
    FLAGS='-DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF -DENABLE_OPENMP=ON' \
    TEST='^(test_config_reference|test_finite_check|test_profiles|test_deterministic_seed|test_orbit_invariants|test_stencil_utils|test_orbit_symplectic_base|test_newton_solver_status|test_spectre_fo_spline_tok2vol|test_spectre_fo_spline_ncsx_freeboundary|test_diag_albert)$'
```

The complete fast suite includes integration cases with 900-second timeouts
(`test/tests/CMakeLists.txt:1087,1098`); it was not run on the shared login node.
The selected tests do not exercise the Python API, MPI, GPU execution or scaling.

Use the QA VMEC file downloaded by CMake's test configuration
(`test/tests/CMakeLists.txt:18–29`) for an eight-particle Boozer run.
The tested `wout.nc` SHA-256 was
`333ab8cef64c4b4f406e76209266fc3bf4d7976e9c28c162983c2120e112e771`:

```bash
mkdir -p run/cineca-smoke
cp examples/hpc/cineca/smoke.in run/cineca-smoke/simple.in
cp test/test_data/wout.nc run/cineca-smoke/wout.nc
/usr/bin/time -p timeout 180 bash -c \
    'cd run/cineca-smoke && ../../build/simple.x simple.in'
```

The input fixes seed 42, `trace_time=1d-5` seconds, `sbeg=0.3`, eight particles
and twenty output steps. Inspect `confined_fraction.dat`, `times_lost.dat` and
`orbit_exit_code.dat` in the run directory. Use a separate run directory for
each concurrent simulation because these output names are fixed.
The observed exit codes were four completed (`0`) and four skipped (`3`):
the default `contr_pp=-1` skips strongly passing particles
(`src/params.f90:60,728`, `src/simple_main.f90:1918–1926`).

For the budget-free batch smoke test, submit from the repository root after
preparing that run directory:

```bash
sbatch examples/hpc/cineca/leonardo-serial.sbatch
squeue -j <job-id>
sacct -j <job-id> --format=JobID,State,ExitCode,Elapsed
```

The template requests one task, four CPUs, 4 GB and three minutes. It sets
`OMP_PLACES=cores`, `OMP_PROC_BIND=close` and `srun --cpu-bind=cores`.
See [Leonardo's partition limits](https://docs.hpc.cineca.it/hpc/leonardo.html#job-managing-and-slurm-partitions).

### Verification record (2026-10-18)

| Check | Result |
| --- | --- |
| Final GNU Release/OpenMP reference rebuild, `cmake --build build --clean-first -j 4` | Passed, 107.24 s wall time |
| Selected `make test-fast` command above | 11/11 passed, 12.99 s CTest time (22.40 s including Make's configure/build) |
| Timed eight-particle example above | Exit 0, 9.29 s wall time, 6,287 field evaluations |
| Output inspection with Python's standard library | 8 particle rows, 8 exit-code rows, 20 confinement rows; all finite; final time `1e-5` s |
| `sbatch examples/hpc/cineca/leonardo-serial.sbatch` | Job 61061586 completed, `sacct` exit `0:0`, 17 s elapsed (9 s application step) |
| GNU/CUDA configure and selected targets below | Passed, 26.91 s configure and 72.28 s compilation; no GPU execution |
| Local worktree GNU 16.2.1 Release/OpenMP | Same 11/11 tests passed in 7.34 s; example exit 0 in 18.98 s |

These are smoke-test timings on shared machines, not performance comparisons.
The final environment script also directs temporary files to `build/tmp`.
One approximately 455 MB CPU reference tree remains at
`$HOME/SIMPLE-efss-cineca/build`; alternate builds were removed. Logs, module
lists, CMake caches and Slurm accounting remain in
`$HOME/SIMPLE-efss-cineca-evidence`. No Fortran or CMake source was changed.

## CUDA and NVHPC

GPU execution requires an allocated GPU. Do not run CUDA tests on login nodes
or the serial partition. The native CUDA path and valid input restrictions are
documented in [DOC/gpu-cuda.md](../../DOC/gpu-cuda.md).

The GNU/CUDA A100 compilation uses the CPU environment above plus CUDA 12.3:

```bash
module load cuda/12.3
cmake -S . -B build-cuda -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DSIMPLE_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DSIMPLE_DETERMINISTIC_FP=ON -DENABLE_PYTHON_INTERFACE=OFF \
    -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF \
    -DLIBNEO_REF=f6090e8150c2f5b1c0afde09baadf79b50b1ffe7
cmake --build build-cuda -j 4 --target \
    simple.x test_cuda_native_rk54.x test_cuda_native_bridge.x
cuobjdump --list-elf build-cuda/libsimple_cuda_native_rk54.so
```

The last command found an `sm_80` cubin. The device code resides in a shared
library (`CMakeLists.txt:120–132`); `ldd build-cuda/simple.x` resolved it and
`libcudart.so.12`. These checks do not launch a GPU kernel.

The NVHPC attempt used this separate build directory and compiler combination
(the explicit temporary directory keeps subsequent builds under `$HOME`):

```bash
module purge
module load profile/base nvhpc/24.5 gcc/12.2.0 cmake/3.27.9 \
    ninja/1.11.1 python/3.11.7 cuda/12.3
export TMPDIR="$HOME/SIMPLE-efss-cineca/build-nvhpc-cuda/tmp"
mkdir -p "$TMPDIR"
cmake -S . -B build-nvhpc-cuda -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_C_COMPILER=nvc \
    -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CUDA_HOST_COMPILER=g++ \
    -DSIMPLE_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DSIMPLE_DETERMINISTIC_FP=OFF -DENABLE_PYTHON_INTERFACE=OFF \
    -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF \
    -DLIBNEO_REF=f6090e8150c2f5b1c0afde09baadf79b50b1ffe7
cmake --build build-nvhpc-cuda -j 4
```

Configuration passed, but compilation failed in Fortnum's
`src/verified/fortnum_interval_qp.f90:22`. Line 10 imports `real128` and line 22
declares `real(qp)`. NVFORTRAN 24.5 reports an illegal negative KIND;
an independent `selected_real_kind(33,4931)` probe returned `-3`.
Fortnum unconditionally includes this file at `src/CMakeLists.txt:181` in the
pinned revision. Upstream main `b29e82cf6e41dc1e5509c2e6722ad99f5a5132bd`
still includes it at line 189. No dependency or physics workaround was applied.
After this fatal error, the remaining compiler process for
`src/generated/fortnum_dop853_tableau.f90` was stopped after more than four
minutes; no successful NVHPC build is claimed.
Also leave deterministic mode off for this attempt: libneo's
`CMakeLists.txt:86–89` explicitly rejects it with non-GNU Fortran.

## Known setup pitfalls

* Missing `netcdf.h` fails at `test/tests/generate_test_boozer_chartmap.c:2`.
  Its target links bare `netcdf` at `test/tests/CMakeLists.txt:466–467`.
  Loading the serial NetCDF-C module supplies the header and library paths.
* Mixing `git/2.38.1` with the OpenMPI stack failed during module loading:
  `/leonardo/prod/opt/modulefiles/base/archive/git/2.38.1:43–44` loads
  `openssh/9.1p1`, which loads the older Kerberos dependency. Use `/usr/bin/git`
  after `module purge`; the minimal serial stack avoids that conflict.
* SIMPLE and libneo main were checked before changing the module recipe:
  their heads were still `fe2de27` and `f6090e8`. Existing compiler warnings
  in fetched Fortplot were nonfatal; this work does not repair them.

## Pitagora: derived recipe, not yet executed

Consult the current [Pitagora guide](https://docs.hpc.cineca.it/hpc/pitagora.html).
The single end-of-task retry, `ssh -o BatchMode=yes -o ConnectTimeout=20
pitagora hostname`, returned status 255 with `Connection closed`. No shell,
module inventory, build or job submission was possible; no access workaround
was attempted. The account-side diagnosis was supplied with the task and was
not independently established by this SSH result.

Use `login-dcgp.pitagora.cineca.it` for AMD Turin CPU builds and
`login-boost.pitagora.cineca.it` for Intel/H100 builds. Inventory modules again;
Leonardo's exact module names are not asserted to exist on Pitagora. Select
GCC, CMake, Ninja, Python, BLAS and NetCDF-C there, save the tested module loads
in an environment script, and export its absolute path as `SIMPLE_ENV`.
Repeat the pinned Release configuration, short test selection and example
above before production use. Rebuild on the destination partition; do not
carry a host-tuned binary across architectures.

For the H100, derive the GNU/CUDA command with
`-DCMAKE_CUDA_ARCHITECTURES=90` in `build-cuda`. This follows NVIDIA's
[Hopper compilation guidance](https://docs.nvidia.com/cuda/hopper-compatibility-guide/index.html#building-applications-with-hopper-architecture-support);
neither compilation nor runtime has been verified on Pitagora.
Prepare `run/cineca-gpu/simple.in` with `gpu_mode='production'`,
`gpu_backend='cuda-native'`, `gpu_method='dopri'`, `gpu_num_devices=1` and
`isw_field_type=2`; keep collisions, walls, classifiers and orbit output disabled.
Copy the equilibrium into that directory. The CPU smoke input does not select
GPU execution.

Partition limits from CINECA, checked 2026-10-18:

| Partition | Maximum nodes/job | Maximum wall time |
| --- | ---: | --- |
| `dcgp_fua_dbg` | 2 | 30 minutes |
| `dcgp_fua_prod` | 64 | 24 hours |
| `boost_fua_dbg` | 2 | 30 minutes |
| `boost_fua_prod` | 16 | 24 hours |

`ptgr_all_serial` is budget-free, CPU-only, limited to four cores and four hours.
The templates [pitagora-cpu.sbatch](../../examples/hpc/cineca/pitagora-cpu.sbatch)
and [pitagora-gpu.sbatch](../../examples/hpc/cineca/pitagora-gpu.sbatch) request
three-minute debug jobs. Replace `PROJECT_ACCOUNT` with an authorized account
before submission; these templates were not submitted in this verification.
The CPU template uses the prepared smoke directory; the GPU template uses
`run/cineca-gpu`. SIMPLE is launched as one process with OpenMP threads.
