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

The build system for SIMPLE is CMake.

Required libraries:
* NetCDF
* LAPACK/BLAS

Supported compilers:
* GNU Fortan
* Intel Fortran

Building:
```bash
cd /path/to/SIMPLE
mkdir build; cd build
cmake ..
make
```
This will produce `libsimple.so` and `simple.x` required to run the code.

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
