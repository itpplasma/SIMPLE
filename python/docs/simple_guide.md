# SIMPLE Python API Guide

Simple functional interface for SIMPLE particle orbit tracing.

## Quick Start

```python
import simple

# Load field once for whole simulation
simple.load_field('wout.nc')

# Sample particles
particles = simple.sample_surface(n_particles=1000, s=0.9)

# Trace orbits
results = simple.trace(particles, tmax=100.0)

# Analyze results
confined = simple.get_confined(results)
lost = simple.get_lost(results)

print(f"Confinement fraction: {confined.shape[1] / particles.shape[1]:.3f}")
```

## Installation

1. **Build SIMPLE with Python support:**
   ```bash
   cd SIMPLE
   make  # Builds pysimple module
   ```

2. **Use simple.py interface:**
   ```python
   import sys
   sys.path.append('/path/to/SIMPLE/python')
   import simple
   ```

## API Reference

### Field Loading

**`simple.load_field(vmec_file)`**
- Loads VMEC equilibrium for whole simulation
- Must be called before sampling or tracing
- Uses existing `init_field()` Fortran function

### Particle Sampling

**`simple.sample_surface(n_particles, s=0.9, **kwargs)`**
- Sample particles on flux surface
- Returns (5, n_particles) array in SoA format
- Uses existing `samplers.f90` functions

**`simple.sample_volume(n_particles, s_inner=0.1, s_outer=0.9, **kwargs)`**
- Sample particles in volume between flux surfaces
- Returns (5, n_particles) array in SoA format

**`simple.load_particles(particle_file, **kwargs)`**
- Load particles from file
- Returns (5, n_particles) array in SoA format

### Orbit Tracing

**`simple.trace(particles, tmax=100.0, integrator='midpoint', **kwargs)`**
- Trace particle orbits using existing Fortran implementation
- `particles`: (5, n_particles) or (n_particles, 5) array
- `tmax`: Maximum integration time
- `integrator`: 'euler', 'midpoint', 'symplectic', 'rk4'
- Returns dictionary with results

### Result Analysis

**`simple.get_confined(results, t_threshold=None)`**
- Get confined particles from simulation results
- Returns (5, n_confined) array

**`simple.get_lost(results, t_threshold=None)`**
- Get lost particles from simulation results
- Returns dict with 'positions' and 'loss_times'

### Utilities

**`simple.info()`**
- Print SIMPLE version and build information

## Example Workflows

### Basic Surface Sampling
```python
import simple

# Setup
simple.load_field('wout.nc')

# Sample and trace
particles = simple.sample_surface(1000, s=0.9)
results = simple.trace(particles, tmax=50.0)

# Analyze
confined = simple.get_confined(results)
print(f"Confined: {confined.shape[1]}/{particles.shape[1]}")
```

### Volume Study
```python
import simple
import numpy as np

simple.load_field('wout.nc')

# Compare different flux surfaces
s_values = np.linspace(0.1, 0.9, 9)
confinement_fractions = []

for s in s_values:
    particles = simple.sample_surface(500, s=s)
    results = simple.trace(particles, tmax=100.0)
    confined = simple.get_confined(results)
    fraction = confined.shape[1] / particles.shape[1]
    confinement_fractions.append(fraction)
    print(f"s={s:.1f}: {fraction:.3f}")
```

### Load from File
```python
import simple

simple.load_field('wout.nc')

# Load particles from file
particles = simple.load_particles('initial_conditions.dat')
results = simple.trace(particles, tmax=200.0, integrator='symplectic')

# Detailed analysis
confined = simple.get_confined(results)
lost = simple.get_lost(results)

print(f"Total: {particles.shape[1]}")
print(f"Confined: {confined.shape[1]}")
print(f"Lost: {lost['positions'].shape[1]}")
print(f"Mean loss time: {np.mean(lost['loss_times']):.2f}")
```

## Notes

- **Field Loading**: `load_field()` must be called before any sampling or tracing
- **SoA Format**: All particle arrays use Structure-of-Arrays (5, n_particles) format
- **Zero-Copy**: Direct access to existing Fortran arrays for performance
- **Thread Safety**: Uses existing OpenMP parallelization from Fortran

## Integration with Existing Code

This API is a thin wrapper around existing `pysimple` module:

```python
# Direct pysimple usage
import pysimple
pysimple.params.vmec_file = 'wout.nc'
pysimple.init_field()
# ... complex setup ...

# Equivalent simple API
import simple
simple.load_field('wout.nc')
particles = simple.sample_surface(1000)
results = simple.trace(particles)
```

The simple API handles all the setup and parameter management automatically.