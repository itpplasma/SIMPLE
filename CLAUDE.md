# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System and Commands

### Building
- **Primary build**: `make` - Uses CMake with Ninja generator to build in `build/` directory
- **Alternative**: `fmp build` - Uses Fortran Package Manager
- **Python package**: Build via scikit-build-core (see pyproject.toml)

### Testing
- **Run tests**: `make test` - Executes ctest in build directory
- Note: Test system is currently under development

### Build Configuration
- Default: Release build (`-DCMAKE_BUILD_TYPE=Release`)
- Available types: Release, Debug, Profile, Fast
- Change config: `make CONFIG=Debug`

### Dependencies
- Required: NetCDF, LAPACK/BLAS, libneo
- Compilers: GNU Fortran or Intel Fortran
- Optional: OpenMP (enabled by default)

## Architecture Overview

SIMPLE is a symplectic particle orbit tracer for fusion plasma physics, designed for stellarator optimization with VMEC equilibria.

### Core Module Structure

**Main Components:**
- `simple_main.f90`: Main simulation orchestration with OpenMP parallelization
- `simple.f90`: Core `Tracer` type and orbit integration interfaces
- `params.f90`: Centralized configuration via namelist input

**Field System** (`src/field/`):
- Abstract `MagneticField` base class in `field_base.f90`
- VMEC implementation in `field_vmec.f90` 
- Canonical coordinates: Boozer, flux, Meiss, Albert variants
- Runtime field selection via function pointers

**Integration System:**
- `orbit_symplectic_base.f90`: Base integrator types and RK coefficients
- `orbit_symplectic.f90`: Multiple symplectic methods (Euler, midpoint, Gauss, Lobatto)
- `orbit_symplectic_quasi.f90`: Quasi-symplectic and RK45 methods

**Coordinate Transformations** (`coordinates/`):
- VMEC ↔ Cylindrical ↔ Cartesian transformations
- Extensible via function pointers

**Particle Management:**
- `samplers.f90`: Surface/volume/file-based particle initialization
- `classification.f90`: Trapped/passing and regular/chaotic orbit classification

### Key Design Patterns

- **Strategy Pattern**: Field evaluation with runtime method selection
- **Template Method**: Common integration interface with specialized implementations
- **Abstract Base**: Extensible field representations

### Data Flow
1. Load VMEC equilibrium and configuration from `simple.in`
2. Initialize particles via selected sampling method
3. Parallel orbit integration (one particle per OpenMP thread)
4. Output confinement statistics and trajectories

### Main Executable
- `simple.x`: Primary executable requiring `simple.in` and VMEC file
- Example inputs in `examples/simple.in` and downloadable `wout.nc`

### Python Interface
- `pysimple` module built automatically with main library
- Located in `python/` directory with integration examples

### Output Files
- `confined_fraction.dat`: Main confinement statistics over time
- `times_lost.dat`: Individual particle loss times and classification
- `start.dat`: Initial conditions (input/output depending on startmode)