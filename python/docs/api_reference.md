# SIMPLE Python API Reference

Complete API documentation for the batch-oriented HPC Python interface to SIMPLE particle orbit tracing.

## Module: simple

### Core Classes

## ParticleBatch

**Class**: `simple.ParticleBatch(n_particles: int)`

Batch container for particle data using existing SoA layout with zero-copy access to `zstart(5, n_particles)` arrays.

### Constructor

```python
ParticleBatch(n_particles: int)
```

**Parameters:**
- `n_particles` (int): Number of particles in the batch. Must be positive.

**Raises:**
- `ValueError`: If `n_particles <= 0`

**Example:**
```python
import simple
particles = simple.ParticleBatch(n_particles=10000)
```

### Properties

#### positions

```python
@property
positions -> np.ndarray
```

Zero-copy view of the underlying `zstart` array with SoA layout.

**Returns:**
- `np.ndarray`: Shape `(5, n_particles)` with:
  - Row 0: s coordinates (flux surface, 0 < s < 1)
  - Row 1: theta coordinates (poloidal angle, 0 to 2π)
  - Row 2: phi coordinates (toroidal angle, 0 to 2π)
  - Row 3: v_par coordinates (parallel velocity)
  - Row 4: mu coordinates (magnetic moment)

**Memory Characteristics:**
- C-contiguous for optimal performance
- Cache-friendly column access pattern
- Zero-copy view of Fortran arrays

**Example:**
```python
particles = simple.ParticleBatch(1000)
pos = particles.positions
print(f"Shape: {pos.shape}")  # (5, 1000)

# Access specific coordinates
s_coords = pos[0, :]      # All s coordinates
first_particle = pos[:, 0]  # All coordinates of first particle
```

#### coordinates

```python
@property
coordinates -> Coordinates
```

Structured access to particle coordinates with named fields.

**Returns:**
- `Coordinates`: Named tuple with fields:
  - `s`: Flux surface coordinates (np.ndarray)
  - `theta`: Poloidal angle coordinates (np.ndarray)
  - `phi`: Toroidal angle coordinates (np.ndarray)
  - `v_par`: Parallel velocity coordinates (np.ndarray)
  - `mu`: Magnetic moment coordinates (np.ndarray)

**Example:**
```python
coords = particles.coordinates
print(f"s range: [{coords.s.min():.3f}, {coords.s.max():.3f}]")
print(f"Mean parallel velocity: {coords.v_par.mean():.3f}")
```

### Initialization Methods

#### initialize_surface

```python
initialize_surface(vmec_file: str, s: float, **kwargs)
```

Initialize particles on a flux surface using deterministic sampling.

**Parameters:**
- `vmec_file` (str): Path to VMEC equilibrium file (.nc format)
- `s` (float): Flux surface coordinate (0 < s < 1)
- `**kwargs`: Additional sampling parameters (reserved for future use)

**Raises:**
- `ValueError`: If s is not in range (0, 1)
- `RuntimeError`: If initialization fails

**Physics:**
- Flux surface s=0 corresponds to magnetic axis
- s=1 corresponds to last closed flux surface
- Particles distributed uniformly in (theta, phi) space
- Velocity and magnetic moment sampled from physical distributions

**Example:**
```python
particles = simple.ParticleBatch(50000)
particles.initialize_surface("wout.nc", s=0.9)

# Verify initialization
coords = particles.coordinates
assert 0.85 < coords.s.min() < 0.95  # s ± 0.05 variation
assert coords.theta.min() >= 0 and coords.theta.max() <= 2*np.pi
```

#### initialize_volume

```python
initialize_volume(vmec_file: str, s_min: float, s_max: float)
```

Initialize particles in a volume between flux surfaces.

**Parameters:**
- `vmec_file` (str): Path to VMEC equilibrium file
- `s_min` (float): Minimum flux surface coordinate
- `s_max` (float): Maximum flux surface coordinate

**Constraints:**
- Must satisfy: 0 < s_min < s_max < 1

**Raises:**
- `ValueError`: If flux surface bounds are invalid
- `RuntimeError`: If initialization fails

**Example:**
```python
particles = simple.ParticleBatch(100000)
particles.initialize_volume("wout.nc", s_min=0.1, s_max=0.95)

coords = particles.coordinates
print(f"s distribution: [{coords.s.min():.3f}, {coords.s.max():.3f}]")
```

#### initialize_from_array

```python
initialize_from_array(particle_data: np.ndarray)
```

Initialize particles from existing coordinate data.

**Parameters:**
- `particle_data` (np.ndarray): Particle coordinate array in one of:
  - SoA format: shape `(5, n_particles)`
  - AoS format: shape `(n_particles, 5)` (automatically transposed)

**Coordinate Order:**
- Index 0: s (flux surface)
- Index 1: theta (poloidal angle)
- Index 2: phi (toroidal angle)
- Index 3: v_par (parallel velocity)
- Index 4: mu (magnetic moment)

**Raises:**
- `ValueError`: If array shape is incompatible

**Example:**
```python
# From SoA format
soa_data = np.random.rand(5, 1000)
soa_data[0, :] = 0.9  # Set all to s=0.9 surface
particles = simple.ParticleBatch(1000)
particles.initialize_from_array(soa_data)

# From AoS format (automatically transposed)
aos_data = np.random.rand(1000, 5)
particles.initialize_from_array(aos_data)
```

### Data Access Methods

#### get_particle

```python
get_particle(index: int) -> np.ndarray
```

Get coordinates for a single particle using cache-friendly column access.

**Parameters:**
- `index` (int): Particle index (0 <= index < n_particles)

**Returns:**
- `np.ndarray`: Shape `(5,)` with `[s, theta, phi, v_par, mu]`

**Performance:**
- Optimized for SoA memory layout
- Cache-friendly column access pattern

**Raises:**
- `IndexError`: If index is out of range

**Example:**
```python
particles = simple.ParticleBatch(1000)
particles.initialize_surface("wout.nc", s=0.9)

# Get first particle coordinates
particle_0 = particles.get_particle(0)
s, theta, phi, v_par, mu = particle_0
print(f"Particle 0: s={s:.3f}, theta={theta:.3f}")
```

#### set_particle

```python
set_particle(index: int, coordinates: np.ndarray)
```

Set coordinates for a single particle.

**Parameters:**
- `index` (int): Particle index
- `coordinates` (np.ndarray): Shape `(5,)` with `[s, theta, phi, v_par, mu]`

**Raises:**
- `IndexError`: If index is out of range
- `ValueError`: If coordinates array has wrong shape

**Example:**
```python
# Modify specific particle
new_coords = np.array([0.85, np.pi, np.pi/2, 0.5, 0.1])
particles.set_particle(42, new_coords)

# Verify modification
modified = particles.get_particle(42)
assert np.allclose(modified, new_coords)
```

### Utility Methods

#### copy

```python
copy() -> ParticleBatch
```

Create a deep copy with independent memory.

**Returns:**
- `ParticleBatch`: Independent copy with same data

**Example:**
```python
original = simple.ParticleBatch(1000)
original.initialize_surface("wout.nc", s=0.9)

# Create independent copy
backup = original.copy()
assert backup.n_particles == original.n_particles
assert not np.shares_memory(backup.positions, original.positions)
```

#### slice

```python
slice(start: int, end: Optional[int] = None) -> ParticleBatch
```

Create a new batch containing a subset of particles.

**Parameters:**
- `start` (int): Starting particle index (inclusive)
- `end` (int, optional): Ending particle index (exclusive). Defaults to `n_particles`.

**Returns:**
- `ParticleBatch`: New batch with copied subset data

**Raises:**
- `ValueError`: If slice range is invalid

**Example:**
```python
particles = simple.ParticleBatch(10000)
particles.initialize_surface("wout.nc", s=0.9)

# Get first 1000 particles
subset = particles.slice(0, 1000)
assert subset.n_particles == 1000

# Get last 500 particles
tail = particles.slice(-500)
assert tail.n_particles == 500
```

#### to_dict / from_dict

```python
to_dict() -> Dict[str, Any]
```

Export data for serialization.

**Returns:**
- `dict`: Serializable representation with keys:
  - `'coordinates'`: Copy of positions array
  - `'n_particles'`: Number of particles
  - `'metadata'`: Backend metadata

```python
@classmethod
from_dict(cls, data: Dict[str, Any]) -> ParticleBatch
```

Reconstruct from serialized data.

**Example:**
```python
# Serialize
particles = simple.ParticleBatch(1000)
particles.initialize_surface("wout.nc", s=0.9)
data = particles.to_dict()

# Reconstruct
restored = simple.ParticleBatch.from_dict(data)
assert np.allclose(restored.positions, particles.positions)
```

---

## BatchResults

**Class**: `simple.BatchResults(backend_results: FortranResultWrapper)`

Results container providing zero-copy access to simulation result arrays with vectorized analysis methods.

### Constructor

Typically created by `trace_orbits()`. For manual construction:

```python
BatchResults(backend_results: FortranResultWrapper)
```

**Parameters:**
- `backend_results`: Wrapper around Fortran result arrays

### Properties

#### loss_times

```python
@property
loss_times -> np.ndarray
```

Access to particle loss times.

**Returns:**
- `np.ndarray`: Shape `(n_particles,)` with loss times
  - Finite values: Time when particle was lost
  - `np.inf`: Particle remained confined

**Physics:**
- Loss occurs when particle exits computational domain
- Typically happens at plasma edge or divertor strike
- Units match simulation time units

**Example:**
```python
results = simple.trace_orbits(particles, tmax=1000.0)
loss_times = results.loss_times

confined_count = np.sum(loss_times == np.inf)
lost_count = np.sum(loss_times < np.inf)
print(f"Confined: {confined_count}, Lost: {lost_count}")
```

#### final_positions

```python
@property
final_positions -> np.ndarray
```

Access to final particle coordinates.

**Returns:**
- `np.ndarray`: Shape `(5, n_particles)` with final coordinates
  - Row 0: final s coordinates
  - Row 1: final theta coordinates
  - Row 2: final phi coordinates
  - Row 3: final v_par coordinates
  - Row 4: final mu coordinates

**Example:**
```python
final_pos = results.final_positions

# Analyze final s distribution
final_s = final_pos[0, :]
print(f"Final s range: [{final_s.min():.3f}, {final_s.max():.3f}]")

# Plot Poincaré section
import matplotlib.pyplot as plt
plt.scatter(final_pos[1, :], final_pos[2, :], s=1, alpha=0.5)
plt.xlabel('θ (poloidal)')
plt.ylabel('φ (toroidal)')
```

#### trap_parameter

```python
@property
trap_parameter -> np.ndarray
```

Access to trapping parameter for orbit classification.

**Returns:**
- `np.ndarray`: Shape `(n_particles,)` with trapping parameters

**Physics:**
- Distinguishes trapped vs passing particles
- trap_parameter > 0: Trapped orbits
- trap_parameter < 0: Passing orbits
- Magnitude indicates orbit characteristics

**Example:**
```python
trap_par = results.trap_parameter
trapped_mask = trap_par > 0
passing_mask = trap_par < 0

print(f"Trapped orbits: {trapped_mask.sum()}")
print(f"Passing orbits: {passing_mask.sum()}")
```

#### perpendicular_invariant

```python
@property
perpendicular_invariant -> np.ndarray
```

Access to perpendicular adiabatic invariant.

**Returns:**
- `np.ndarray`: Shape `(n_particles,)` with perpendicular invariants

**Physics:**
- Conserved quantity in guiding center motion
- Related to magnetic moment and field strength
- Conservation indicates quality of symplectic integration

#### confined_mask / lost_mask

```python
@property
confined_mask -> np.ndarray
@property 
lost_mask -> np.ndarray
```

Boolean masks for particle classification.

**Returns:**
- `np.ndarray`: Shape `(n_particles,)` boolean arrays
  - `confined_mask`: True for confined particles (`loss_time == np.inf`)
  - `lost_mask`: True for lost particles (`loss_time < np.inf`)

**Example:**
```python
confined = results.confined_mask
lost = results.lost_mask

# Subset operations
confined_positions = results.final_positions[:, confined]
lost_loss_times = results.loss_times[lost]
```

### Analysis Methods

#### confinement_statistics

```python
confinement_statistics(time_bins: int = 50) -> ConfinementStats
```

Comprehensive vectorized confinement analysis.

**Parameters:**
- `time_bins` (int): Number of histogram bins for loss time distribution

**Returns:**
- `ConfinementStats`: Data class with fields:
  - `n_total` (int): Total number of particles
  - `n_confined` (int): Number of confined particles
  - `n_lost` (int): Number of lost particles
  - `confined_fraction` (float): Fraction remaining confined (0 to 1)
  - `mean_loss_time` (float): Mean loss time for lost particles
  - `loss_time_distribution` (Tuple[np.ndarray, np.ndarray]): Histogram data

**Performance:**
- Vectorized NumPy operations
- O(n_particles) complexity
- Memory-efficient for large datasets

**Example:**
```python
stats = results.confinement_statistics(time_bins=100)

print(f"Confinement Summary:")
print(f"  Total: {stats.n_total:,} particles")
print(f"  Confined: {stats.n_confined:,} ({stats.confined_fraction:.1%})")
print(f"  Lost: {stats.n_lost:,}")
if stats.n_lost > 0:
    print(f"  Mean loss time: {stats.mean_loss_time:.3f}")

# Plot loss time distribution
hist_counts, hist_bins = stats.loss_time_distribution
import matplotlib.pyplot as plt
plt.bar(hist_bins[:-1], hist_counts, width=np.diff(hist_bins))
plt.xlabel('Loss Time')
plt.ylabel('Number of Particles')
```

#### get_confined_particles / get_lost_particles

```python
get_confined_particles() -> np.ndarray
get_lost_particles() -> Tuple[np.ndarray, np.ndarray, np.ndarray]
```

Extract subsets for confined/lost particles.

**Returns:**
- `get_confined_particles()`: Final positions `(5, n_confined)`
- `get_lost_particles()`: Tuple of:
  - Final positions `(5, n_lost)`
  - Loss times `(n_lost,)`
  - Trap parameters `(n_lost,)`

**Example:**
```python
# Analyze confined particles
confined_final = results.get_confined_particles()
if confined_final.shape[1] > 0:
    confined_s = confined_final[0, :]
    print(f"Confined particles: s ∈ [{confined_s.min():.3f}, {confined_s.max():.3f}]")

# Analyze lost particles
lost_final, lost_times, lost_trap = results.get_lost_particles()
if len(lost_times) > 0:
    early_loss = lost_times < 100.0  # Lost early
    print(f"Early losses: {early_loss.sum()} particles")
```

#### analyze_by_surface

```python
analyze_by_surface(s_bins: np.ndarray) -> Dict[str, np.ndarray]
```

Analyze confinement statistics by flux surface.

**Parameters:**
- `s_bins` (np.ndarray): Flux surface bin edges

**Returns:**
- `dict`: Analysis results with keys:
  - `'s_centers'`: Bin center coordinates
  - `'confined_fraction'`: Fraction confined in each bin
  - `'mean_loss_time'`: Mean loss time in each bin
  - `'particle_count'`: Number of particles in each bin

**Example:**
```python
# Create 20 flux surface bins
s_bins = np.linspace(0, 1, 21)
surface_analysis = results.analyze_by_surface(s_bins)

# Print confinement profile
for s, frac, count in zip(surface_analysis['s_centers'],
                         surface_analysis['confined_fraction'],
                         surface_analysis['particle_count']):
    if count > 0:
        print(f"s={s:.2f}: {frac:.1%} confined ({count:,} particles)")
```

### Export Methods

#### save_hdf5

```python
save_hdf5(filename: str, append: bool = False)
```

Export to HDF5 format with compression.

**Parameters:**
- `filename` (str): Output HDF5 file path
- `append` (bool): If True, append to existing file

**Features:**
- Gzip compression for efficiency
- Metadata attributes included
- Supports multiple batch appending

**Requires:**
- `h5py` package installed

**Example:**
```python
# Single export
results.save_hdf5("simulation_results.h5")

# Multi-run accumulation
for run in range(5):
    # ... run simulation ...
    results.save_hdf5("multi_run.h5", append=True)
```

#### save_numpy / load_numpy

```python
save_numpy(filename: str)

@classmethod
load_numpy(cls, filename: str) -> BatchResults
```

NumPy compressed format export/import.

**Example:**
```python
# Export
results.save_numpy("results.npz")

# Import
loaded = simple.BatchResults.load_numpy("results.npz")
assert loaded.n_particles == results.n_particles
```

#### summary

```python
summary() -> str
```

Generate human-readable summary.

**Returns:**
- `str`: Multi-line summary with key statistics

**Example:**
```python
print(results.summary())
# Output:
# BatchResults Summary:
#   Total particles: 100,000
#   Confined: 85,432 (85.4%)
#   Lost: 14,568 (14.6%)
#   Mean loss time: 234.567
```

---

## Simulation Functions

### trace_orbits

```python
trace_orbits(
    particles: ParticleBatch,
    tmax: float,
    integrator: str = 'symplectic_midpoint',
    openmp_threads: Optional[int] = None,
    memory_efficient: bool = False,
    config: Optional[Dict[str, Any]] = None,
    verbose: bool = False
) -> BatchResults
```

Main batch orbit tracing function leveraging existing OpenMP parallelization.

**Performance Guarantees:**
- <5% overhead vs direct Fortran execution
- Zero-copy access to existing SoA arrays
- Leverages existing thread-per-particle OpenMP pattern

**Parameters:**
- `particles` (ParticleBatch): Initialized particle batch
- `tmax` (float): Maximum simulation time
- `integrator` (str): Integration method name
- `openmp_threads` (int, optional): OpenMP thread count (None = auto)
- `memory_efficient` (bool): Enable memory optimization for large datasets
- `config` (dict, optional): Pre-created configuration dictionary
- `verbose` (bool): Print detailed execution information

**Available Integrators:**
- `'symplectic_euler'`: First-order symplectic
- `'symplectic_midpoint'`: Second-order symplectic (recommended)
- `'symplectic_gauss'`: Higher-order Gauss method
- `'symplectic_lobatto'`: Lobatto symplectic method
- `'quasi_symplectic'`: Quasi-symplectic comparison
- `'rk45'`: Adaptive Runge-Kutta (non-symplectic)

**Returns:**
- `BatchResults`: Simulation results with zero-copy array access

**Example:**
```python
# Basic usage
particles = simple.ParticleBatch(10000)
particles.initialize_surface("wout.nc", s=0.9)

results = simple.trace_orbits(
    particles,
    tmax=1000.0,
    integrator='symplectic_midpoint',
    openmp_threads=8,
    verbose=True
)

# Advanced configuration
config = simple.create_configuration(
    vmec_file="wout.nc",
    tmax=1000.0,
    integrator='symplectic_gauss',
    dtau=0.1
)

results = simple.trace_orbits(
    particles,
    config=config,
    memory_efficient=True
)
```

### create_configuration

```python
create_configuration(
    vmec_file: str,
    tmax: float,
    integrator: str = 'symplectic_midpoint',
    dtau: Optional[float] = None,
    field_type: str = 'vmec',
    **kwargs
) -> Dict[str, Any]
```

Create simulation configuration dictionary with automatic parameter validation.

**Parameters:**
- `vmec_file` (str): Path to VMEC equilibrium file
- `tmax` (float): Maximum simulation time
- `integrator` (str): Integration method name
- `dtau` (float, optional): Time step (auto-calculated if None)
- `field_type` (str): Magnetic field type
- `**kwargs`: Additional simulation parameters

**Time Step Auto-Calculation:**
- `symplectic_euler`: tmax/10000 (conservative)
- `symplectic_midpoint`: tmax/5000 (moderate)
- `symplectic_gauss/lobatto`: tmax/2000 (larger for higher-order)
- `rk45`: tmax/1000 (adaptive stepping)

**Returns:**
- `dict`: Configuration for `trace_orbits()`

**Raises:**
- `ValueError`: If integrator name is unknown
- `FileNotFoundError`: If VMEC file doesn't exist

**Example:**
```python
# Automatic time step
config = simple.create_configuration(
    vmec_file="wout.nc",
    tmax=1000.0,
    integrator='symplectic_midpoint'
)

# Manual time step
config = simple.create_configuration(
    vmec_file="wout.nc",
    tmax=1000.0,
    integrator='symplectic_gauss',
    dtau=0.05
)

# Use configuration
results = simple.trace_orbits(particles, config=config)
```

### Convenience Functions

#### quick_simulation

```python
quick_simulation(
    n_particles: int,
    vmec_file: str = "wout.nc",
    s_surface: float = 0.9,
    tmax: float = 1000.0,
    integrator: str = 'symplectic_midpoint'
) -> BatchResults
```

One-line simulation with sensible defaults.

**Example:**
```python
# Rapid prototyping
results = simple.quick_simulation(
    n_particles=10000,
    vmec_file="wout.nc",
    s_surface=0.9,
    tmax=1000.0
)
print(results.summary())
```

#### parameter_sweep

```python
parameter_sweep(
    base_config: Dict[str, Any],
    parameter_name: str,
    parameter_values: list,
    particles: ParticleBatch,
    parallel: bool = False
) -> Dict[Any, BatchResults]
```

Perform parameter sweep with batch processing.

**Example:**
```python
base_config = simple.create_configuration(
    vmec_file="wout.nc",
    tmax=1000.0
)

# Sweep integration time step
dtau_values = [0.05, 0.1, 0.2, 0.5]
sweep_results = simple.parameter_sweep(
    base_config,
    'dtau',
    dtau_values,
    particles
)

# Analyze convergence
for dtau, result in sweep_results.items():
    stats = result.confinement_statistics()
    print(f"dtau={dtau}: {stats.confined_fraction:.4f}")
```

---

## Large-Scale Processing

### process_large_simulation

```python
process_large_simulation(
    vmec_file: str,
    n_total: int,
    tmax: float,
    batch_size: int = 100_000,
    output_file: str = 'results.h5',
    s_surface: float = 0.9,
    integrator: str = 'symplectic_midpoint',
    verbose: bool = True,
    memory_limit_gb: Optional[float] = None
) -> StreamResults
```

Process millions of particles with constant memory usage through batch streaming.

**Memory Usage:** Constant - independent of total particle count  
**Performance:** Linear scaling with existing OpenMP optimization  
**Output:** Streaming HDF5 for large dataset handling

**Parameters:**
- `vmec_file` (str): Path to VMEC equilibrium file
- `n_total` (int): Total number of particles to simulate
- `tmax` (float): Maximum simulation time
- `batch_size` (int): Particles per batch (controls memory usage)
- `output_file` (str): Output HDF5 file path
- `s_surface` (float): Flux surface for particle initialization
- `integrator` (str): Integration method
- `verbose` (bool): Print progress information
- `memory_limit_gb` (float, optional): Memory limit in GB

**Returns:**
- `StreamResults`: Handle for accessing streaming results

**Example:**
```python
# Process 10 million particles
stream_results = simple.process_large_simulation(
    vmec_file="wout.nc",
    n_total=10_000_000,
    tmax=1000.0,
    batch_size=100_000,
    output_file='large_sim.h5',
    s_surface=0.9,
    verbose=True,
    memory_limit_gb=16.0
)

# Get overall statistics
summary = stream_results.get_summary_statistics()
print(f"Total confined fraction: {summary['confined_fraction']:.3f}")
```

### StreamResults

**Class**: `simple.StreamResults(output_file: str, total_particles: int)`

Container for streaming simulation results with efficient access to large datasets.

#### get_summary_statistics

```python
get_summary_statistics() -> Dict[str, Any]
```

Get summary statistics without loading all data into memory.

**Returns:**
- `dict`: Summary statistics with keys:
  - `'total_particles'`: Total number of particles
  - `'total_confined'`: Number of confined particles  
  - `'total_lost'`: Number of lost particles
  - `'confined_fraction'`: Overall confinement fraction
  - `'mean_loss_time'`: Mean loss time across all batches
  - `'median_loss_time'`: Median loss time

**Example:**
```python
summary = stream_results.get_summary_statistics()
print(f"Processed {summary['total_particles']:,} particles")
print(f"Confined fraction: {summary['confined_fraction']:.3f}")
```

#### load_batch_results / iter_batch_results

```python
load_batch_results(batch_index: int) -> BatchResults
iter_batch_results() -> Iterator[BatchResults]
```

Access individual batch results or iterate over all batches.

**Example:**
```python
# Load specific batch
batch_0 = stream_results.load_batch_results(0)
stats_0 = batch_0.confinement_statistics()

# Iterate over all batches
for i, batch_results in enumerate(stream_results.iter_batch_results()):
    stats = batch_results.confinement_statistics()
    print(f"Batch {i}: {stats.confined_fraction:.3f} confined")
```

#### combine_all_results

```python
combine_all_results() -> BatchResults
```

Combine all batch results into single `BatchResults` object.

**Warning:** Loads all data into memory - may require significant RAM.

**Example:**
```python
# Only for datasets that fit in memory
combined = stream_results.combine_all_results()
overall_stats = combined.confinement_statistics()
```

---

## Memory Management

### Memory Utilities

#### estimate_memory_usage

```python
estimate_memory_usage(
    n_particles: int,
    include_results: bool = True,
    overhead_factor: float = 1.5
) -> Dict[str, float]
```

Estimate memory usage for given particle count.

**Example:**
```python
memory_1M = simple.estimate_memory_usage(1_000_000)
print(f"1M particles: {memory_1M['total_gb']:.1f} GB")
print(f"Breakdown: {memory_1M['particles_mb']:.0f} MB particles, "
      f"{memory_1M['results_mb']:.0f} MB results")
```

#### optimize_batch_size

```python
optimize_batch_size(
    target_memory_gb: float,
    n_total: int,
    safety_factor: float = 0.8
) -> int
```

Optimize batch size for memory constraints.

**Example:**
```python
optimal_batch = simple.optimize_batch_size(
    target_memory_gb=8.0,
    n_total=5_000_000
)
print(f"Optimal batch size: {optimal_batch:,} particles")
```

### MemoryMonitor

```python
simple.MemoryMonitor(name: str = "simulation", verbose: bool = True)
```

Context manager for monitoring memory usage.

**Example:**
```python
with simple.MemoryMonitor("large_simulation") as monitor:
    results = simple.trace_orbits(particles, tmax=1000.0)
    current_memory = monitor.sample()
    
stats = monitor.get_stats()
print(f"Memory increase: {stats['increase_gb']:.1f} GB")
```

---

## Performance and Validation

### Performance Benchmarking

#### benchmark_performance

```python
benchmark_performance(
    n_particles: int,
    tmax: float = 100.0,
    integrator: str = 'symplectic_midpoint',
    n_runs: int = 3,
    openmp_threads: Optional[int] = None
) -> Dict[str, float]
```

Benchmark Python API performance.

**Returns:**
- `dict`: Performance metrics:
  - `'mean_time'`: Mean execution time
  - `'particles_per_second'`: Processing rate
  - `'overhead_estimate'`: Estimated API overhead

**Example:**
```python
metrics = simple.benchmark_performance(
    n_particles=100_000,
    tmax=100.0,
    n_runs=5
)

print(f"Performance: {metrics['particles_per_second']:.0f} particles/sec")
print(f"API overhead: {metrics['overhead_estimate']:.1%}")
```

### Golden Record Validation

#### validate_golden_record

```python
validate_golden_record(
    reference_results: Union[str, BatchResults],
    test_particles: ParticleBatch,
    config: Dict[str, Any],
    rtol: float = 1e-10,
    atol: float = 1e-12
) -> Dict[str, bool]
```

Validate results against reference implementation.

**Example:**
```python
# Load reference results
reference = simple.BatchResults.load_numpy("reference.npz")

# Run test
validation = simple.validate_golden_record(
    reference,
    test_particles,
    config,
    rtol=1e-10
)

for test, passed in validation.items():
    print(f"{test}: {'PASS' if passed else 'FAIL'}")
```

### Results Comparison

#### compare_results

```python
simple.compare_results(
    results1: BatchResults,
    results2: BatchResults,
    rtol: float = 1e-10,
    atol: float = 1e-12
) -> Dict[str, bool]
```

Compare two `BatchResults` objects element-wise.

**Example:**
```python
comparison = simple.compare_results(results1, results2)
if comparison['all_match']:
    print("Results match within tolerance")
else:
    for array, matches in comparison.items():
        if not matches:
            print(f"Mismatch in {array}")
```

---

## Data Classes

### ConfinementStats

**Data Class**: Comprehensive confinement analysis results.

**Fields:**
- `n_total` (int): Total number of particles
- `n_confined` (int): Number of confined particles
- `n_lost` (int): Number of lost particles
- `confined_fraction` (float): Fraction confined (0 to 1)
- `mean_loss_time` (float): Mean loss time for lost particles
- `loss_time_distribution` (Tuple): Histogram `(counts, bins)`

### Coordinates

**Data Class**: Structured particle coordinate access.

**Fields:**
- `s` (np.ndarray): Flux surface coordinates
- `theta` (np.ndarray): Poloidal angle coordinates
- `phi` (np.ndarray): Toroidal angle coordinates
- `v_par` (np.ndarray): Parallel velocity coordinates
- `mu` (np.ndarray): Magnetic moment coordinates

---

## Error Handling

### Common Exceptions

- `ValueError`: Invalid parameter values, array shapes, coordinate ranges
- `RuntimeError`: Backend initialization failures, simulation execution errors
- `FileNotFoundError`: Missing VMEC files or HDF5 datasets
- `IndexError`: Particle index out of range
- `ImportError`: Missing optional dependencies (h5py, psutil)

### Example Error Handling

```python
try:
    particles = simple.ParticleBatch(1000000)
    particles.initialize_surface("missing_file.nc", s=0.9)
    results = simple.trace_orbits(particles, tmax=1000.0)
    
except FileNotFoundError as e:
    print(f"VMEC file not found: {e}")
except ValueError as e:
    print(f"Invalid parameter: {e}")
except RuntimeError as e:
    print(f"Simulation failed: {e}")
except ImportError as e:
    print(f"Missing dependency: {e}")
```

---

## Integration with Scientific Python

### NumPy Integration

All arrays are NumPy arrays supporting standard operations:

```python
# Statistical analysis
final_s = results.final_positions[0, :]
s_mean = np.mean(final_s)
s_std = np.std(final_s)

# Boolean indexing
confined_mask = results.confined_mask
confined_positions = results.final_positions[:, confined_mask]

# Vectorized operations
normalized_loss_times = results.loss_times / results.loss_times.max()
```

### Matplotlib Visualization

```python
import matplotlib.pyplot as plt

# Confinement statistics
stats = results.confinement_statistics()
hist_counts, hist_bins = stats.loss_time_distribution

plt.figure(figsize=(12, 4))

# Loss time histogram
plt.subplot(1, 2, 1)
plt.bar(hist_bins[:-1], hist_counts, width=np.diff(hist_bins))
plt.xlabel('Loss Time')
plt.ylabel('Number of Particles')
plt.title('Loss Time Distribution')

# Poincaré plot
plt.subplot(1, 2, 2)
final_pos = results.final_positions
confined = results.confined_mask
lost = results.lost_mask

plt.scatter(final_pos[1, confined], final_pos[2, confined], 
           s=1, alpha=0.6, label='Confined', color='blue')
plt.scatter(final_pos[1, lost], final_pos[2, lost], 
           s=1, alpha=0.6, label='Lost', color='red')
plt.xlabel('θ (poloidal)')
plt.ylabel('φ (toroidal)')
plt.legend()
plt.title('Poincaré Section')

plt.tight_layout()
plt.show()
```

### HDF5 Integration

```python
import h5py

# Manual HDF5 analysis
with h5py.File('results.h5', 'r') as f:
    # List all batch groups
    batch_groups = [name for name in f.keys() if name.startswith('batch_')]
    
    for group_name in batch_groups:
        grp = f[group_name]
        n_particles = grp.attrs['n_particles']
        confined_frac = grp.attrs['confined_fraction']
        print(f"{group_name}: {confined_frac:.3f} confined ({n_particles:,} particles)")
```

This completes the comprehensive API reference for the SIMPLE Python batch-oriented HPC interface. All functions include type hints, detailed parameters, examples, and physics context appropriate for fusion plasma research.