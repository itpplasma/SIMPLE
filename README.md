[![DOI](https://zenodo.org/badge/146476544.svg)](https://zenodo.org/badge/latestdoi/146476544)

# SIMPLE
**S**ymplectic **I**ntegration **M**ethods for **P**article **L**oss **E**stimation

SIMPLE computes statistical losses of guiding-center orbits for particles of given mass, charge and
energy from the volume of 3D magnetic configurations. Orbits are traced via a symplectic integrator [1,2]
that guarantees conservation of invariants of motion within fixed bounds over long integration periods.
A classifier based on Poincarè plots [1] allows for accelerated prediction for confinement of regular (non-chaotic) orbits.

The main focus of SIMPLE is to make computation of fusion alpha particle losses fast enough to be used directly in
stellarator optimization codes. For that reason currently 3D configurations with nested magnetic flux surfaces are
supported based on a VMEC equilibrium file in NetCDF format, and orbits are computed without taking collisions into account.

The code is free to use and modify under the MIT License and links to Runge-Kutta-Fehlberg routines in
`SRC/contrib/rkf45.f90` from https://people.sc.fsu.edu/~jburkardt/f_src/rkf45/rkf45.html under the GNU LGPL License.

## Building
Run `make` to produce a `build` directory including the main executable
`simple.x`, main library `libsimple.so` and Python module `pysimple`.

Required libraries:
* NetCDF
* LAPACK/BLAS

Required compiler:
* GNU Fortran (gfortran)

For Python wrappers, do
```bash
pip install git+https://github.com/jameskermode/f90wrap
pip install -e . --no-build-isolation
```

### Python API

SIMPLE provides a clean module-level Python API for orbit tracing and classification.

#### Basic Usage

```python
import pysimple

# Initialize with VMEC file and parameters
pysimple.init('wout.nc', deterministic=True, trace_time=1e-3)

# Sample particles on a flux surface
particles = pysimple.sample_surface(100, s=0.5)

# Trace orbits in parallel
results = pysimple.trace_parallel(particles)
print(f"Lost: {(results['loss_times'] < 1e-3).sum()} particles")
```

#### Initialization

**`pysimple.init(vmec_file, **params)`**

Initialize SIMPLE with a VMEC equilibrium file and optional parameters.

Parameters:
- `vmec_file`: Path to VMEC NetCDF file (wout.nc)
- `deterministic=True`: Use fixed random seed for reproducibility
- `trace_time=1e-3`: Integration time in seconds
- `ntestpart=100`: Number of test particles
- `npoiper2=64`: Integration steps per poloidal transit
- `integmode`: Integration method (default: MIDPOINT)
- `isw_field_type`: Field type (0=TEST, 2=VMEC, 3=BOOZER, etc.)
- Any other Fortran parameter from params.f90

#### Particle Sampling

**`pysimple.sample_surface(n_particles, s)`**

Sample particles uniformly on a flux surface.

Returns: `(5, n_particles)` array with columns `[s, theta, phi, p_abs, v_par]`

**`pysimple.sample_volume(n_particles, s_inner, s_outer)`**

Sample particles uniformly in a volume between two flux surfaces.

**`pysimple.load_particles(filename)`**

Load particles from a text file (e.g., start.dat format).

#### Orbit Tracing

**`pysimple.trace_parallel(positions, integrator='midpoint')`**

Trace multiple particle orbits in parallel.

Parameters:
- `positions`: `(5, n_particles)` array of initial conditions
- `integrator`: Integration method ('midpoint', 'rk45', 'gauss2', etc.)

Returns dictionary:
- `'final_positions'`: `(5, n_particles)` final positions
- `'loss_times'`: `(n_particles,)` loss times
- `'trap_parameter'`: `(n_particles,)` trapping parameter
- `'perpendicular_invariant'`: `(n_particles,)` perpendicular invariant

**`pysimple.trace_orbit(position, integrator='midpoint', return_trajectory=False)`**

Trace a single particle orbit.

If `return_trajectory=True`, returns full trajectory arrays.

#### Orbit Classification

**`pysimple.classify_parallel(positions, integrator='midpoint')`**

Classify particle orbits (trapped/passing, regular/chaotic).

Returns dictionary including all trace_parallel outputs plus:
- `'passing'`: `(n_particles,)` boolean (True=passing, False=trapped)
- `'lost'`: `(n_particles,)` boolean
- `'jpar'`: `(n_particles,)` J-parallel conservation (0-2)
- `'topology'`: `(n_particles,)` topological classification (0-2)
- `'minkowski'`: `(n_particles,)` Minkowski dimension (0=unclassified, 1=regular, 2=chaotic)

Example:
```python
pysimple.init('wout.nc', tcut=0.1, deterministic=True, trace_time=1e-3)
particles = pysimple.sample_surface(1000, s=0.5)
classified = pysimple.classify_parallel(particles)

regular = classified['minkowski'] == 1
chaotic = classified['minkowski'] == 2
trapped = ~classified['passing']

print(f"Regular: {regular.sum()}, Chaotic: {chaotic.sum()}")
print(f"Trapped: {trapped.sum()}, Passing: {(~trapped).sum()}")
```

#### Integrator Constants

Available integrators (from `orbit_symplectic_base.f90`):
- `pysimple.RK45` (0): Runge-Kutta 4/5
- `pysimple.MIDPOINT` (3): Symplectic midpoint (default)
- `pysimple.GAUSS1` (4): Gauss 1st order
- `pysimple.GAUSS2` (5): Gauss 2nd order
- `pysimple.LOBATTO3` (15): Lobatto 3rd order

#### Complete Examples

See `examples/simple_api.py` for complete working examples.

## Usage

SIMPLE currently runs on a single node with OpenMP shared memory parallelization with one particle per thread and background
fields residing in main memory. The main executable is `simple.x`.

### Input
* `simple.in`
* a VMEC NetCDF equlibrium (wout.nc) file with name specified in `simple.in`

An example input file with explanation of each parameter can be found in `examples/simple.in`. An example `wout.nc` can be obtained by
```bash
wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc
```

In addition `start.dat` is either an input for given (`startmode>=2`) or an output (`startmode` 0 or 1) for randomly generated initial conditions.
Diagnostics for slow convergence of Newton iterations are written in `fort.6601`.

### Output
`confined_fraction.dat` is the main output, containing four columns:
1. Physical time
2. Confined fraction of passing particles
3. Confined fraction of trapped particles
4. Total number of particles

The sum of 2. and 3. yields the overall confined fraction at each time.

`times_lost.dat` contains the loss time of each particle. Columns are:
1. Particle index. Corresponds to line number in start.dat .
2. Time t_loss [s] when the particle is lost. Possible values are: -1, `trace_time`, or any other value between 0 and trace_time.
  * If never lost or classified as regular, maximum tracing time `trace_time` is written.
  * If ignored due to contr_pp, which defines deep passing region as confined (we don consider them anymore), -1 is written.
3. Trapping parameter trap_par that is 1 for deeply trapped, 0 for tp boundary
   and negative for passing. Eq. (3.1) in Accelerated Methods paper.
   Whenever trap_par < contr_pp, particle is not traced and counted as confined.

### Comparing Commits ("Golden Record")
To execute test_against_legacy_behaviour.py, execute `get_test_data.sh`, run SIMPLE and create the following symlinks in the `tests` directory after running SIMPLE in the `test_data` directory:
* ../../../test_data/wout.nc --> wout.nc
*  ../../../../examples/simple_full.in --> simple.in
*  ../../../../build/start.dat start.dat

## References
When using this code for scientific publications, please cite the according references:

[1] C. G. Albert, S. V. Kasilov, and W. Kernbichler,
Accelerated methods for direct computation of fusion alpha particle losses within stellarator optimization. J. Plasma Phys 86, 815860201 (2020), https://doi.org/10.1017/S0022377820000203

[2] C. G. Albert, S. V. Kasilov, and W. Kernbichler,
Symplectic integration with non-canonical quadrature for guiding-center orbits in magnetic confinement devices. J. Comp. Phys 403, 109065 (2020), https://doi.org/10.1016/j.jcp.2019.109065, preprint on https://arxiv.org/abs/1903.06885
