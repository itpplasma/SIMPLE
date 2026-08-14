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
modes are documented in `DOC/config.md` and `examples/simple_full.in`. The
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

For a Python example, see `examples/simple_api.py`:

```python
import pysimple

pysimple.init("wout.nc", deterministic=True, trace_time=1e-3)
particles = pysimple.sample_surface(100, s=0.5)
results = pysimple.trace_parallel(particles)
lost = ((results["loss_times"] > 0.0) &
        (results["loss_times"] < 1e-3)).sum()
skipped = (results["loss_times"] < 0.0).sum()
print(f"Lost: {lost} particles (skipped deep-passing: {skipped})")
```

Examples for classification, trajectories, and plotting are also included. The
Python interface accepts the parameters documented in `src/params.f90` and
`DOC/config.md`.

## CPU quickstart

Copy `examples/simple.in` and provide a VMEC file named `wout.nc`:

```bash
mkdir -p run
cp examples/simple.in run/simple.in
wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/\
wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O run/wout.nc
(cd run && ../build/simple.x simple.in)
```

`examples/simple_full.in` contains the complete input file with comments.
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

Numerically unresolved orbits are reported separately from physical losses.

## GPU builds and runs

GPU runs are configured in `simple.in`. No GPU environment variables are
needed. Set `gpu_mode = 'production'` and use Boozer guiding-center tracing
(`isw_field_type = 2`) without collisions, walls, orbit output, or classifiers:

```fortran
gpu_mode = 'production'
gpu_backend = 'cuda-native'
gpu_method = 'dopri'       ! tuned Dormand-Prince
gpu_landreman = .True.     ! segmented alpha-loss objective [3]
```

The native CUDA path requires the CUDA toolkit, a C++17 compiler, and a normal
GNU Fortran build:

```bash
cmake -S . -B build-cuda -G Ninja \
  -DSIMPLE_ENABLE_CUDA=ON
cmake --build build-cuda -j
```

With the `simple.in` settings above, run:

```bash
(cd run && ../build-cuda/simple.x simple.in)
```

Native CUDA supports tuned DOPRI and Cash-Karp; `gpu_method = 'dopri'` is the
default. The tolerance is `relerr` from the same `simple.in` file. Set
`gpu_start_coordinates = 'boozer'` only when `start.dat` already contains
Boozer coordinates. Set `gpu_particle_profile` to write the per-particle CSV.
The Landreman controls are `gpu_t_block`, `gpu_loss_tau`, `gpu_maxloss`, and
`gpu_min_timestep`; see [3] for the related loss-threshold objective.

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
