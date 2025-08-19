# SIMPLE Python Interface

Simple functional API for SIMPLE particle orbit tracer.

## Usage

```python
import simple

# Load field
simple.load_field('wout.nc')

# Sample particles  
particles = simple.sample_surface(1000, s=0.3)

# Trace orbits
results = simple.trace(particles, tmax=0.1, integrator=simple.EXPL_IMPL_EULER)

# Analyze
confined = simple.get_confined(results)
print(f"Confined: {confined.shape[1]}/{particles.shape[1]}")
```

## Functions

- `load_field(vmec_file)` - Load VMEC equilibrium
- `sample_surface(n_particles, s=0.3)` - Sample on flux surface
- `sample_volume(n_particles, s_inner=0.1, s_outer=0.9)` - Sample in volume
- `load_particles(file)` - Load from file
- `trace(particles, tmax=0.1, integrator=EXPL_IMPL_EULER)` - Trace orbits
- `get_confined(results)` - Get confined particles
- `get_lost(results)` - Get lost particles

## Integration Methods

- `simple.RK45 = 0`
- `simple.EXPL_IMPL_EULER = 1` (default)
- `simple.IMPL_EXPL_EULER = 2`
- `simple.MIDPOINT = 3`
- `simple.GAUSS1 = 4`, `simple.GAUSS2 = 5`, `simple.GAUSS3 = 6`, `simple.GAUSS4 = 7`
- `simple.LOBATTO3 = 15`

See `example_simple_usage.py` for complete examples.