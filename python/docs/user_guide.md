# SIMPLE Python API User Guide

High-performance Python interface for SIMPLE particle orbit tracing with batch-oriented HPC capabilities and zero-copy access to existing Fortran SoA data structures.

## Quick Start

### Basic Particle Batch Processing

```python
import simple

# Create particle batch
particles = simple.ParticleBatch(n_particles=10000)

# Initialize particles on flux surface
particles.initialize_surface("wout.nc", s=0.9)

# Run simulation
results = simple.trace_orbits(particles, tmax=1000.0)

# Analyze results
stats = results.confinement_statistics()
print(f"Confined fraction: {stats.confined_fraction:.3f}")
```

## Installation

### Prerequisites

1. **Build SIMPLE with Python support:**
   ```bash
   cd SIMPLE
   make  # Builds with Python support automatically
   ```

2. **Verify pysimple availability:**
   ```python
   import sys
   sys.path.append('build')
   import pysimple  # Should work without errors
   ```

### Installing the Batch API

The batch-oriented API is available immediately after building SIMPLE:

```python
# Add to your Python path
import sys
sys.path.append('/path/to/SIMPLE/python')

import simple
print(simple.__version__)  # Should print "1.0.0"
```

### VMEC Equilibrium File

Download a test VMEC file for examples:

```bash
wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc
```

## Core Concepts

### Batch-Oriented Architecture

The API is designed for HPC workflows processing thousands to millions of particles:

- **ParticleBatch**: Zero-copy wrapper around existing `zstart(5, n_particles)` arrays
- **BatchResults**: Performance-optimized access to simulation results  
- **trace_orbits()**: Leverages existing OpenMP parallelization with <5% overhead

### Structure of Arrays (SoA) Layout

SIMPLE uses Structure of Arrays for optimal memory performance:

```python
# Particle positions: shape (5, n_particles)
# positions[0, :] = s coordinates (flux surface)
# positions[1, :] = theta coordinates (poloidal angle)  
# positions[2, :] = phi coordinates (toroidal angle)
# positions[3, :] = v_par coordinates (parallel velocity)
# positions[4, :] = mu coordinates (magnetic moment)
```

This layout provides:
- Cache-friendly column access: `positions[:, particle_id]`
- Contiguous memory for GPU transfer
- Zero-copy access to existing Fortran arrays

## Particle Initialization

### Surface Sampling

Initialize particles on a flux surface:

```python
particles = simple.ParticleBatch(n_particles=50000)
particles.initialize_surface(
    vmec_file="wout.nc", 
    s=0.9  # Flux surface coordinate (0 < s < 1)
)

# Access coordinates directly
coords = particles.coordinates
print(f"s range: [{coords.s.min():.3f}, {coords.s.max():.3f}]")
```

### Volume Sampling

Initialize particles in a volume:

```python
particles = simple.ParticleBatch(n_particles=100000)
particles.initialize_volume(
    vmec_file="wout.nc",
    s_min=0.1,
    s_max=0.95
)
```

### Custom Initialization

Initialize from existing arrays:

```python
# From SoA format
particle_data = np.random.rand(5, 1000)
particles = simple.ParticleBatch(1000)
particles.initialize_from_array(particle_data)

# From AoS format (automatically transposed)
aos_data = np.random.rand(1000, 5)
particles.initialize_from_array(aos_data)
```

## Simulation Execution

### Basic Simulation

```python
results = simple.trace_orbits(
    particles,
    tmax=1000.0,
    integrator='symplectic_midpoint',
    openmp_threads=8,
    verbose=True
)
```

### Available Integrators

- `'symplectic_euler'`: First-order symplectic (conservative time step)
- `'symplectic_midpoint'`: Second-order symplectic (recommended default)
- `'symplectic_gauss'`: Higher-order Gauss method
- `'symplectic_lobatto'`: Lobatto symplectic method
- `'quasi_symplectic'`: Quasi-symplectic for comparison
- `'rk45'`: Adaptive Runge-Kutta (non-symplectic)

### Advanced Configuration

```python
config = simple.create_configuration(
    vmec_file="wout.nc",
    tmax=1000.0,
    integrator='symplectic_gauss',
    dtau=0.1,  # Time step (auto-calculated if None)
    field_type='vmec'
)

results = simple.trace_orbits(
    particles,
    config=config,
    memory_efficient=True
)
```

## Results Analysis

### Confinement Statistics

```python
stats = results.confinement_statistics()
print(f"Total particles: {stats.n_total}")
print(f"Confined: {stats.n_confined} ({stats.confined_fraction:.1%})")
print(f"Lost: {stats.n_lost}")
print(f"Mean loss time: {stats.mean_loss_time:.3f}")

# Loss time distribution
hist_counts, hist_bins = stats.loss_time_distribution
```

### Accessing Result Arrays

Zero-copy access to all result data:

```python
# Loss times (np.inf for confined particles)
loss_times = results.loss_times

# Final particle positions (5, n_particles)
final_positions = results.final_positions

# Physics quantities
trap_parameters = results.trap_parameter
perp_invariants = results.perpendicular_invariant

# Boolean masks
confined_mask = results.confined_mask
lost_mask = results.lost_mask
```

### Particle Subsets

```python
# Get confined particles only
confined_final = results.get_confined_particles()

# Get lost particle data
lost_final, lost_times, lost_trap = results.get_lost_particles()
```

### Surface Analysis

Analyze confinement by flux surface:

```python
s_bins = np.linspace(0, 1, 21)  # 20 flux surface bins
surface_analysis = results.analyze_by_surface(s_bins)

print("Confinement by surface:")
for i, (s, frac) in enumerate(zip(surface_analysis['s_centers'], 
                                 surface_analysis['confined_fraction'])):
    print(f"  s={s:.2f}: {frac:.1%} confined")
```

## Large-Scale Simulations

### Memory-Efficient Streaming

For millions of particles, use streaming to maintain constant memory usage:

```python
# Process 10 million particles in batches of 100k
stream_results = simple.process_large_simulation(
    vmec_file="wout.nc",
    n_total=10_000_000,
    tmax=1000.0,
    batch_size=100_000,
    output_file='results.h5',
    s_surface=0.9,
    verbose=True
)

# Get summary without loading all data
summary = stream_results.get_summary_statistics()
print(f"Overall confined fraction: {summary['confined_fraction']:.3f}")
```

### Memory Usage Estimation

```python
# Estimate memory for different particle counts
memory_1M = simple.estimate_memory_usage(1_000_000)
print(f"1M particles: {memory_1M['total_gb']:.1f} GB")

# Optimize batch size for available memory
optimal_batch = simple.optimize_batch_size(
    target_memory_gb=8.0,
    n_total=5_000_000
)
print(f"Optimal batch size: {optimal_batch:,} particles")
```

### Streaming Results Analysis

```python
# Load specific batch results
batch_0 = stream_results.load_batch_results(0)
batch_stats = batch_0.confinement_statistics()

# Iterate over all batches
for i, batch_results in enumerate(stream_results.iter_batch_results()):
    stats = batch_results.confinement_statistics()
    print(f"Batch {i}: {stats.confined_fraction:.3f} confined")

# Combine all results (requires memory for full dataset)
combined = stream_results.combine_all_results()
```

## Performance Optimization

### OpenMP Threading

```python
# Control OpenMP threads explicitly
results = simple.trace_orbits(
    particles,
    tmax=1000.0,
    openmp_threads=16  # Match your CPU cores
)

# Or set environment variable
import os
os.environ['OMP_NUM_THREADS'] = '16'
```

### Memory Monitoring

```python
# Monitor memory usage during simulation
with simple.MemoryMonitor("large_simulation") as monitor:
    results = simple.trace_orbits(particles, tmax=1000.0)
    print(f"Peak memory: {monitor.peak_memory:.1f} GB")
```

### Performance Benchmarking

```python
# Benchmark the Python API performance
metrics = simple.benchmark_performance(
    n_particles=100_000,
    tmax=100.0,
    n_runs=3
)

print(f"Mean simulation time: {metrics['mean_time']:.3f}s")
print(f"Particles per second: {metrics['particles_per_second']:.0f}")
print(f"Estimated overhead: {metrics['overhead_estimate']:.1%}")
```

## Data Export and Analysis

### HDF5 Export

```python
# Export to HDF5 for large datasets
results.save_hdf5("simulation_results.h5")

# Append multiple runs
for run in range(5):
    # ... run simulation ...
    results.save_hdf5("multi_run.h5", append=True)
```

### NumPy Export

```python
# Export to compressed NumPy format
results.save_numpy("results.npz")

# Load back for analysis
loaded_results = simple.BatchResults.load_numpy("results.npz")
```

### Integration with Scientific Python

```python
import matplotlib.pyplot as plt
import numpy as np

# Plot loss time distribution
stats = results.confinement_statistics()
hist_counts, hist_bins = stats.loss_time_distribution

plt.figure(figsize=(10, 6))
plt.bar(hist_bins[:-1], hist_counts, width=np.diff(hist_bins), alpha=0.7)
plt.xlabel('Loss Time')
plt.ylabel('Number of Particles')
plt.title('Particle Loss Time Distribution')
plt.show()

# Poincaré plot of final positions
final_pos = results.final_positions
confined_mask = results.confined_mask

plt.figure(figsize=(12, 5))

# Confined particles
plt.subplot(1, 2, 1)
plt.scatter(final_pos[1, confined_mask], final_pos[2, confined_mask], 
           s=1, alpha=0.5, label='Confined')
plt.xlabel('θ (poloidal angle)')
plt.ylabel('φ (toroidal angle)')
plt.title('Confined Particles')

# Lost particles
plt.subplot(1, 2, 2)
lost_mask = results.lost_mask
plt.scatter(final_pos[1, lost_mask], final_pos[2, lost_mask], 
           s=1, alpha=0.5, color='red', label='Lost')
plt.xlabel('θ (poloidal angle)')
plt.ylabel('φ (toroidal angle)')
plt.title('Lost Particles')

plt.tight_layout()
plt.show()
```

## Parameter Sweeps

```python
# Define base configuration
base_config = simple.create_configuration(
    vmec_file="wout.nc",
    tmax=1000.0,
    integrator='symplectic_midpoint'
)

# Sweep time step
dtau_values = [0.05, 0.1, 0.2, 0.5]
sweep_results = simple.parameter_sweep(
    base_config,
    parameter_name='dtau',
    parameter_values=dtau_values,
    particles=particles
)

# Analyze convergence
for dtau, result in sweep_results.items():
    stats = result.confinement_statistics()
    print(f"dtau={dtau}: {stats.confined_fraction:.4f} confined")
```

## Quick Simulations

For rapid testing and prototyping:

```python
# One-line simulation with sensible defaults
results = simple.quick_simulation(
    n_particles=10000,
    vmec_file="wout.nc",
    s_surface=0.9,
    tmax=1000.0
)
```

## Integration with Existing Workflow

### Compatibility with Fortran SIMPLE

The Python API uses the same underlying implementation:

```python
# Same results as running: ./simple.x with equivalent simple.in
python_results = simple.trace_orbits(particles, tmax=1000.0)
```

### Golden Record Validation

Validate against reference implementations:

```python
# Compare with reference results
reference_results = simple.BatchResults.load_numpy("reference.npz")
validation = simple.validate_golden_record(
    reference_results,
    test_particles,
    config,
    rtol=1e-10
)

for test, passed in validation.items():
    print(f"{test}: {'PASS' if passed else 'FAIL'}")
```

## Troubleshooting

### Common Issues

**"pysimple module not available"**
- Ensure SIMPLE was built with `make` (includes Python support)
- Check that `build/` directory contains `pysimple.*` files
- Verify Python path includes SIMPLE build directory

**Memory errors with large simulations**
- Use `process_large_simulation()` for millions of particles
- Optimize batch size with `optimize_batch_size()`
- Monitor memory usage with `MemoryMonitor`

**Performance issues**
- Set `OMP_NUM_THREADS` to match CPU cores
- Use SoA-optimized access patterns (column access over row access)
- Enable `memory_efficient=True` for large datasets

**Array shape warnings**
- Normal during initialization - arrays are automatically resized
- Indicates zero-copy integration with existing Fortran implementation

### Performance Validation

```python
# Validate SoA performance characteristics
from simple.core.batch import validate_soa_performance

performance = validate_soa_performance(particles, verbose=True)
if not performance['cache_friendly']:
    print("Warning: SoA layout may not be optimal")
```

### Debugging Integration

```python
# Check backend integration status
backend = simple.backends.fortran.get_backend()
print(f"Backend initialized: {backend._initialized}")

# Inspect array metadata
arrays = particles._arrays
metadata = arrays.get_metadata()
print(f"Available arrays: {[arr['name'] for arr in metadata['available_arrays']]}")
```

## Next Steps

- **API Reference**: See `api_reference.md` for complete function documentation
- **Examples**: Check `python/examples/` for working code samples
- **Performance**: Review benchmarking results and optimization strategies
- **Advanced Features**: Explore streaming, parameter sweeps, and HDF5 integration

The SIMPLE Python API provides production-ready HPC capabilities while maintaining the scientific accuracy and performance of the underlying Fortran implementation.