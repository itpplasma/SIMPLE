# SIMPLE Python Interface

Simple functional API for SIMPLE particle orbit tracer.

## Usage

```python
import simple

# Load field
simple.load_field('wout.nc')

# Sample particles  
particles = simple.sample_surface(1000, s=0.9)

# Trace orbits
results = simple.trace(particles, tmax=100.0, integrator=simple.MIDPOINT)

# Analyze
confined = simple.get_confined(results)
print(f"Confined: {confined.shape[1]}/{particles.shape[1]}")
```

## Functions

- `load_field(vmec_file)` - Load VMEC equilibrium
- `sample_surface(n_particles, s=0.9)` - Sample on flux surface
- `sample_volume(n_particles, s_inner=0.1, s_outer=0.9)` - Sample in volume
- `load_particles(file)` - Load from file
- `trace(particles, tmax=100.0, integrator=MIDPOINT)` - Trace orbits
- `get_confined(results)` - Get confined particles
- `get_lost(results)` - Get lost particles

## Constants

- `simple.EULER = 1`
- `simple.MIDPOINT = 2` 
- `simple.SYMPLECTIC = 3`
- `simple.RK4 = 4`

See `example_simple_usage.py` for complete examples.