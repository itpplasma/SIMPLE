[![DOI](https://zenodo.org/badge/146476544.svg)](https://zenodo.org/badge/latestdoi/146476544)

# SIMPLE
**S**ymplectic **I**ntegration **M**ethods for **P**article **L**oss **E**stimation

SIMPLE computes guiding-center orbit losses for particles in three-dimensional
magnetic configurations. It traces orbits with symplectic integrators and can
classify regular and chaotic orbits from Poincare plots. The main application
is fast evaluation of fusion alpha-particle losses during stellarator
optimization.

The standard run uses a VMEC equilibrium in NetCDF format, nested magnetic
flux surfaces, and collisionless guiding-center orbits. Other field and orbit
modes are documented in `examples/simple_full.in` and `DOC/config.md`. The
software is distributed under the MIT License. The FIRM3D-compatible tracing
implementation has separate provenance information in
[LICENSES/FIRM3D-MIT.txt](LICENSES/FIRM3D-MIT.txt).
The segmented loss objective follows the related GPU tracing and optimization
protocol described by Landreman et al. [3].

## Building

The default build produces `build/simple.x`, the SIMPLE library, and, when
the Python prerequisites are available, the `pysimple` module:

```bash
make
```

Build requirements:

- CMake 3.22 or newer, with Ninja.
- GNU Fortran and Git.
- NetCDF-C and NetCDF-Fortran.
- LAPACK and BLAS.
- OpenMP.

On common systems:

```bash
# Arch Linux
sudo pacman -S gcc-fortran cmake ninja netcdf netcdf-fortran lapack blas

# Ubuntu or Debian
sudo apt install gfortran cmake ninja-build libnetcdf-dev libnetcdff-dev \
    liblapack-dev libblas-dev

# macOS with Homebrew
brew install gcc cmake ninja netcdf netcdf-fortran lapack
```

To select a libneo branch, tag, or commit, pass it explicitly:

```bash
make LIBNEO_REF=my-branch
```

To use a local libneo checkout:

```bash
make LIBNEO_PATH=/path/to/libneo
```

For a reproducible floating-point build and the fast test suite:

```bash
make build-deterministic
make test-fast
```

### Python interface

Create the recommended repository-local environment with:

```bash
./setup-venv.sh
source .venv/bin/activate
```

The Python package `pysimple` runs the same code as `simple.x` and is
configured with the same `simple.in`:

```python
import pysimple

pysimple.init(config="simple.in")            # or pysimple.init("wout.nc", trace_time=1e-3)
result = pysimple.trace_parallel(pysimple.sample())
```

See `python/README.rst` and `examples/simple_api.py`.

### Changed-line lint

Pull requests check only added/modified lines, so legacy formatting outside the
change does not block review. Install the two local tools and run the same check from
the repository root (replace `origin/main` with another base revision when
needed):

```bash
python -m pip install fprettify==0.3.7 ruff==0.14.1
python tools/check_changed_lines.py --base origin/main
```

The checker compares tracked working files (including staged and unstaged edits)
against the merge base; stage new files with `git add` first. Free-form Fortran is
checked against `fprettify` output using `.fprettify`. Python files under
`python/` and `tools/` use Ruff lint and format checks. Formatter edits and lint
diagnostics are filtered to changed source lines; files are never rewritten.
Run the checker's regression tests with:

```bash
python -m unittest discover -s test/python -p test_changed_line_lint.py -v
```

CI keeps the existing GNU test jobs and adds GCC 13 smoke tests plus GCC 14
fast-suite coverage on non-draft PRs and pushes to `main`. Coverage reports are
attached as the `fast-suite-coverage` artifact and summarized in the job summary.
To collect coverage locally with GNU Fortran, GCC's matching `gcov`, and `lcov`:

```bash
make test-fast CONFIG=Debug \
    FLAGS='-DENABLE_COVERAGE=ON -DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=OFF'
cmake --build build --target coverage
lcov --summary build/coverage_filtered.info
```

The `CI` workflow's **Run workflow** (`workflow_dispatch`) also offers
non-blocking Intel ifx and LLVM flang smoke jobs. These reuse the Make/CMake
build with Python integration and GNU-only deterministic options disabled;
compiler portability failures remain visible for investigation. The additional
GCC jobs run on manual dispatch only when `full` is `true`, preserving the small
libneo smoke gate.

### EUFUS mirror

On pushes to `main` and tags in `itpplasma/SIMPLE`, `mirror-eufus.yml` pushes
`main` and all tags to `https://gitlab.eufus.psnc.pl/itpplasma/simple.git` using
`EUFUS_GITLAB_TOKEN` with username `oauth2`. With no token, it succeeds with a
notice. Forks and pull requests cannot run the mirror. Pushes use ordinary
fast-forward/tag checks; divergent mirror refs require maintainer attention.

## CPU quickstart

Copy `examples/simple.in` and provide a VMEC file named `wout.nc`:

```bash
mkdir -p run
cp examples/simple.in run/simple.in
wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/\
wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O run/wout.nc
(cd run && ../build/simple.x simple.in)
```

`examples/simple_full.in` lists every input key with its default value and
meaning; a run needs only the keys that differ.
`examples/run_example.sh` performs the same setup in `/tmp/simple_example`.

The main outputs are:

- `confined_fraction.dat`: time, passing fraction, trapped fraction, and
  resolved-particle count.
- `unresolved_fraction.dat`: time, unresolved fraction, and total particle
  count.
- `times_lost.dat`: particle index, loss time, trapping parameter, initial
  radial coordinate, perpendicular invariant, and final state.
- `orbit_exit_code.dat`: particle index, exit code, loss time, and boundary
  diagnostics.
- `results.nc`: per-orbit results and diagnostics when NetCDF output is
  enabled.

The trapping parameter (`trap_par` in `times_lost.dat` and `results.nc`,
`trap_parameter` in `pysimple`) is evaluated at the start point as

    trap_par = (mu*B_max/E - 1) * B_min/(B_max - B_min)
             = ((1 - lambda**2)*B_max/B - 1) * B_min/(B_max - B_min),

with pitch `lambda = v_par/v` and `B_min`, `B_max` on the start flux surface.
Negative values are passing, 0 is the trapped-passing boundary and 1 is a
deeply trapped particle reflecting at `B_min`. The passing/trapped split in
`confined_fraction.dat` uses `trap_par < 0`; the fractions there are divided
by the resolved-particle count at each time.

Numerically unresolved orbits are reported separately from physical losses.

## GPU quickstart

The supported production GPU backend is native CUDA. Configure it in
`simple.in` for Boozer guiding-center tracing (`isw_field_type = 2`), without
collisions, walls, orbit output, or classifiers:

```fortran
gpu_mode = 'production'
gpu_backend = 'cuda-native'
gpu_method = 'dopri'       ! tuned Dormand-Prince
gpu_landreman = .True.     ! segmented alpha-loss objective [3]
```

Build with the CUDA toolkit, CMake, Ninja, and either GNU or NVIDIA Fortran.
Set the CUDA architecture for the target GPU. Perlmutter's A100 uses `80`:

```bash
cmake -S . -B build-cuda -G Ninja \
  -DSIMPLE_ENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=80
cmake --build build-cuda -j
```

With the `simple.in` settings above, run:

```bash
(cd run && ../build-cuda/simple.x simple.in)
```

The complete GPU configuration and controls are in
[DOC/gpu-cuda.md](DOC/gpu-cuda.md). Cluster-specific build and Slurm examples
for Perlmutter, aCluster, and sCluster are in [RUN/README.md](RUN/README.md).

## Boozer chartmaps

SIMPLE can read precomputed Boozer chartmaps from GVEC, booz_xform, and other
field sources:

```bash
python tools/gvec_to_boozer_chartmap.py parameter_final.ini State_final.dat \
    boozer_chartmap.nc
python tools/booz_xform_to_boozer_chartmap.py boozmn_case.nc \
    boozer_chartmap.nc
```

Use `examples/simple_chartmap.in` as the input template. Set both `field_input`
and `coord_input` to the chartmap. No VMEC file is required at runtime.
`sbeg` is normalized toroidal flux, while the chartmap radial coordinate is
`rho = sqrt(s)`. The exporter handles the documented Boozer sign convention.
See `docs/boozer-chartmap-schema.rst` for the NetCDF schema and scaling rules.

## Golden-record tests

To compare numerical output with a reference revision:

```bash
make test-golden-main
```

The golden-record scripts build the selected reference and compare its output
with the current build. Use `make test-all` for the complete local suite.

## References

[1] C. G. Albert, S. V. Kasilov, and W. Kernbichler, "Accelerated methods for
direct computation of fusion alpha particle losses within stellarator
optimization," *Journal of Plasma Physics* 86, 815860201 (2020),
https://doi.org/10.1017/S0022377820000203

[2] C. G. Albert, S. V. Kasilov, and W. Kernbichler, "Symplectic integration
with non-canonical quadrature for guiding-center orbits in magnetic confinement
devices," *Journal of Computational Physics* 403, 109065 (2020),
https://doi.org/10.1016/j.jcp.2019.109065

[3] M. Landreman, M. Czekanski, A. Giuliani, B. Jang, and R. Conlin,
"Bayesian optimization of stellarator alpha-particle confinement using
data-informed parameter spaces and dimensionality reduction," arXiv:2606.19523
(2026), https://arxiv.org/abs/2606.19523

When using the FIRM3D-compatible tracing path, also retain the provenance in
[LICENSES/FIRM3D-MIT.txt](LICENSES/FIRM3D-MIT.txt).
