# SIMPLE
**S**ymplectic **I**ntegration **M**ethods for **P**article **L**oss **E**stimation

SIMPLE computes statistical losses of guiding-center orbits for particles of given mass, charge and 
energy from the volume of 3D magnetic configurations. Orbits are traced via a symplectic integrator [1] 
that guarantees conservation of invariants of motion within fixed bounds over long integration periods.
A classifier based on Poincarè plots allows for accelerated prediction for confinement of regular (non-chaotic) orbits.

The main focus of SIMPLE is to make computation of fusion alpha particle losses fast enough to be used directly in 
stellarator optimization codes. For that reason currently 3D configurations with nested magnetic flux surfaces are 
supported based on a VMEC equilibrium file in NetCDF format, and orbits are computed without taking collisions into account. 

## Building

## Usage

SIMPLE currently runs on a single node with OpenMP shared memory parallelization with one particle per thread and background
fields residing in main memory. As an input it takes
* `simple.in`
* a VMEC NetCDF equlibrium (wout.nc) file with name specified in `simple.in`

An example input file with explanation of each parameter can be found in `examples/simple.in`.

The main output is `confined_fraction.dat`, containing four columns:
1) Physical time
2) Confined fraction of passing particles
3) Confined fraction of trapped particles
4) Total number of particles

The sum of 2) and 3) yields the overall confined fraction at each time.

In addition `start.dat` is either an input for given or an output for randomly generated initial conditions. Diagnostics for slow convergence of Newton
iterations are written in `fort.6601`.

## References
When using this code for scientific publications, please cite

[1] C. G. Albert, S. V. Kasilov, and W. Kernbichler, 
Symplectic integration with non-canonical quadrature for guiding-center orbits in magnetic confinement devices. J. Comp. Phys (forthcoming, 2019), arXiv:1903.06885
