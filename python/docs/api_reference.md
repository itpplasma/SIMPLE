# SIMPLE Python API Reference

Pure interface layer providing zero-copy access to existing SIMPLE Fortran functionality.

## Design Principle

The Python API is a **pure interface layer** that:
- Delegates all computation to existing, proven Fortran implementations
- Provides zero-copy access to existing Fortran arrays
- Adds <5% overhead vs direct Fortran execution
- Implements no independent computational logic

## Core Modules

### simple.ParticleBatch

Zero-copy wrapper around existing `zstart(5, n_particles)` arrays.

```python
batch = simple.ParticleBatch(n_particles=1000)
batch.initialize_from_samplers("wout.nc", method='surface')

# Zero-copy access to Fortran arrays
positions = batch.positions  # Shape (5, n_particles)
coords = batch.coordinates   # Structured access: s, theta, phi, v_par, mu
```

**Methods:**
- `initialize_from_samplers(vmec_file, method, **kwargs)` - Use existing samplers.f90
- `initialize_from_array(data)` - Initialize from array data
- `copy()` - Create independent copy
- `slice(start, end)` - Create subset batch

### simple.trace_orbits()

Direct interface to existing `simple_main.f90` execution.

```python
results = simple.trace_orbits(
    particles,
    tmax=500.0,
    integrator='symplectic_midpoint',
    openmp_threads=4
)
```

**Parameters:**
- `particles`: ParticleBatch with initialized data
- `tmax`: Simulation time
- `integrator`: Method name (maps to existing Fortran integrators)
- `openmp_threads`: Thread count (sets OMP_NUM_THREADS)
- `config`: Pre-created configuration dict

### simple.BatchResults

Zero-copy wrapper around existing Fortran result arrays.

```python
stats = results.confinement_statistics()
loss_times = results.loss_times  # Direct access to times_lost array
final_pos = results.final_positions  # Direct access to zend array
```

**Properties:**
- `loss_times` - Zero-copy view of existing `times_lost` array
- `final_positions` - Zero-copy view of existing `zend` array
- `confined_mask` - Boolean mask for confined particles
- `lost_mask` - Boolean mask for lost particles

## Samplers Interface

Direct access to existing `samplers.f90` functions.

### simple.SurfaceSampler

```python
sampler = simple.SurfaceSampler("wout.nc")
particle_data = sampler.sample_surface_fieldline(n_particles=1000)
```

**Methods:**
- `sample_surface_fieldline(n_particles)` → `samplers.f90::sample_surface_fieldline()`
- `sample_grid(grid_density)` → `samplers.f90::sample_grid()`
- `sample_volume_single(n, s_inner, s_outer)` → `samplers.f90::sample_volume_single()`

### simple.VolumeSampler

```python
sampler = simple.VolumeSampler("wout.nc")
particle_data = sampler.sample_volume_single(1000, s_inner=0.1, s_outer=0.9)
```

### simple.FileSampler

```python
sampler = simple.FileSampler("wout.nc")
particle_data = sampler.load_from_file("start.dat")  # → sample_read()
ants_data = sampler.sample_points_ants()            # → sample_points_ants()
```

## Memory Utilities

Large-scale processing coordination using existing batch capabilities.

### simple.process_large_simulation()

```python
results = simple.process_large_simulation(
    total_particles=1_000_000,
    vmec_file="wout.nc",
    tmax=500.0,
    max_memory_mb=2000.0
)
```

Coordinates multiple calls to existing `trace_orbits()` without implementing independent streaming logic.

### Utility Functions

- `estimate_memory_usage(n_particles)` - Based on existing Fortran array sizes
- `optimize_batch_size(target_memory_mb)` - Calculate optimal batch size
- `run_million_particle_simulation()` - Convenience for large-scale processing

## Configuration

### Integrator Mapping

Python names map directly to existing Fortran `integmode` values:

```python
_INTEGRATOR_MAP = {
    'symplectic_euler': 1,
    'symplectic_midpoint': 2,
    'symplectic_gauss': 3,
    'symplectic_lobatto': 4,
    'quasi_symplectic': 5,
    'rk45': 6
}
```

### Configuration Creation

```python
config = simple.create_configuration(
    vmec_file="wout.nc",
    tmax=500.0,
    integrator='symplectic_midpoint',
    dtau=1e-6  # Optional - uses Fortran defaults if None
)
```

## Convenience Functions

### Quick Simulation

```python
results = simple.quick_simulation(
    n_particles=1000,
    vmec_file="wout.nc",
    s_surface=0.9,
    tmax=500.0
)
```

### Batch Creation

```python
# Surface sampling
batch = simple.create_surface_batch(1000, "wout.nc", s=0.9)

# Volume sampling  
batch = simple.create_volume_batch(1000, "wout.nc", s_inner=0.1, s_outer=0.9)

# File loading
batch = simple.load_batch_from_file("wout.nc", "start.dat")
```

## Performance Validation

### Benchmark Interface Overhead

```python
metrics = simple.benchmark_performance(
    n_particles=10000,
    tmax=100.0,
    integrator='symplectic_midpoint',
    n_runs=3
)
print(f"Particles/sec: {metrics['particles_per_second']}")
```

### Golden Record Validation

```python
comparison = simple.validate_golden_record(
    reference_file="reference_results.npz",
    test_particles=particles,
    config=config
)
```

## Array Memory Layout

All arrays use existing Fortran SoA (Structure of Arrays) layout:

```
zstart(5, n_particles):
  zstart(1, :) = s coordinates (flux surface)
  zstart(2, :) = theta coordinates (poloidal angle)  
  zstart(3, :) = phi coordinates (toroidal angle)
  zstart(4, :) = v_par coordinates (parallel velocity)
  zstart(5, :) = mu coordinates (magnetic moment)
```

Zero-copy access ensures cache-friendly column access patterns optimized for HPC workflows.

## Error Handling

All validation delegated to existing Fortran implementations:
- Parameter validation through existing namelist processing
- Physical constraints through existing field evaluation
- Numerical stability through existing integrator implementations

Python layer provides only basic argument checking and file existence validation.