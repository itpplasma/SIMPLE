# SIMPLE KISS Interface

**Simple functions, not OOP classes - Direct access to existing Fortran functionality**

This interface addresses user feedback about over-engineered Java-style OOP complexity by providing a minimal functional API that directly calls proven Fortran implementations.

## Key Principles

- ✅ **KISS Design**: Simple functions instead of complex classes
- ✅ **Direct Fortran Integration**: Real calls to `samplers.f90` and `simple_main.f90`
- ✅ **Global VMEC Loading**: Load field once like `simple.x`, not per-sampler
- ✅ **Zero Reimplementation**: Pure wrapper around existing proven code
- ✅ **Minimal API**: Expose essential functionality without abstraction layers

## Quick Start

```python
import simple_kiss as simple

# Simple functional interface - exactly as requested
particles = simple.sample_surface('wout.nc', n_particles=1000, s=0.9)
results = simple.trace_particles('wout.nc', particles, tmax=100.0)
confined = simple.get_confined(results)

print(f"Confinement fraction: {simple.get_confinement_fraction(results):.2%}")
```

## Core Functions

### Particle Sampling (Direct calls to `samplers.f90`)

```python
# Surface sampling using sample_surface_fieldline()
particles = simple.sample_surface('wout.nc', n_particles=1000, s=0.9)

# Volume sampling using sample_volume_single() 
particles = simple.sample_volume('wout.nc', n_particles=500, s_inner=0.1, s_outer=0.9)

# File loading using sample_read()
particles = simple.load_particles('wout.nc', 'start.dat')
```

### Orbit Tracing (Direct calls to `simple_main.f90`)

```python
# Full orbit simulation using main simulation loop
results = simple.trace_particles('wout.nc', particles, tmax=100.0, 
                                integmode=2, dtau=1e-6)

# Results contain:
# - 'times_lost': Loss times (inf = confined)
# - 'final_coords': Final particle coordinates  
# - 'trap_par': Trapped/passing classification
# - 'perp_inv': Perpendicular adiabatic invariant
```

### Analysis Utilities

```python
# Simple analysis functions
confined = simple.get_confined(results)              # Boolean array
fraction = simple.get_confinement_fraction(results)  # Float 0-1

# One-function complete simulation
results = simple.quick_simulation('wout.nc', n_particles=1000, tmax=100.0)
```

## Architecture

### VMEC Loading Strategy
- **Problem**: Old API loaded VMEC per-sampler (inefficient, wrong architecture)  
- **Solution**: Load VMEC once globally like `simple.x` workflow
- **Implementation**: `_ensure_vmec_loaded()` manages global state

### Fortran Integration
- **Real Integration**: Actual calls to `pysimple.samplers` and `pysimple.simple_main`
- **Zero-Copy Arrays**: Direct views of Fortran SoA arrays
- **Proven Code**: No reimplementation - uses existing tested functions

### Comparison with Old API

| Old (Over-engineered) | New (KISS) |
|---------------------|------------|
| `SurfaceSampler(vmec).sample_surface_fieldline(n)` | `sample_surface(vmec, n)` |
| `Simulation(tracer, batch).run()` | `trace_particles(vmec, particles, tmax)` |
| Multiple classes, unclear Fortran integration | Simple functions, direct Fortran calls |
| VMEC loaded per-sampler | VMEC loaded once globally |

## Backend Requirements

The KISS interface requires the `pysimple` module built via f90wrap:

```bash
# Build from project root
make

# Check backend availability
python -c "import simple_kiss; print(simple_kiss.get_backend_info())"
```

## Examples

- `examples/kiss_demo.py` - Complete demonstration
- `examples/basic_usage.py` - Integration with existing examples

## Design Philosophy

This interface follows the KISS (Keep It Simple, Stupid) principle by:

1. **Eliminating unnecessary abstraction** - Direct function calls instead of class hierarchies
2. **Matching Fortran workflow** - Same VMEC loading and simulation structure as `simple.x`
3. **Minimal API surface** - Only essential functions, no feature creep
4. **Zero reimplementation** - Pure wrapper around proven Fortran code

The result is a clean, simple interface that directly exposes the existing high-performance Fortran implementations without OOP complexity.