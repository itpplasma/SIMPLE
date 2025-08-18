# SIMPLE Python API Examples

Comprehensive examples demonstrating the batch-oriented HPC Python interface for SIMPLE particle orbit tracing.

## Overview

These examples showcase the key capabilities of the SIMPLE Python API:

- **Batch-oriented processing**: Zero-copy access to existing Fortran SoA arrays
- **HPC performance**: <5% overhead vs direct Fortran execution  
- **Memory efficiency**: Constant memory usage for millions of particles
- **Scientific analysis**: Physics-focused tools for fusion plasma research

## Examples

### 1. [basic_batch_processing.py](basic_batch_processing.py)

**Purpose**: Fundamental usage patterns and API validation

**Key Features:**
- Particle batch creation and initialization
- Surface and volume sampling
- Basic simulation execution
- SoA performance validation
- Integrator comparison
- Data access patterns
- Export capabilities

**Run Command:**
```bash
python basic_batch_processing.py
```

**Prerequisites:**
- SIMPLE built with Python support (`make`)
- Internet connection (downloads test VMEC file)

**Expected Output:**
- Performance metrics and timing
- Confinement statistics
- SoA memory layout validation
- Integrator comparison results

---

### 2. [large_scale_streaming.py](large_scale_streaming.py)

**Purpose**: Memory-efficient processing of millions of particles

**Key Features:**
- Streaming simulation with constant memory usage
- Memory estimation and batch size optimization
- HDF5 output for large datasets
- Memory monitoring and profiling
- Performance scaling analysis

**Run Command:**
```bash
python large_scale_streaming.py
```

**Prerequisites:**
- h5py: `pip install h5py`
- psutil: `pip install psutil`
- At least 4GB RAM recommended

**Expected Output:**
- Memory usage optimization
- Streaming simulation progress
- Batch-by-batch analysis
- Performance scaling metrics

---

### 3. [performance_comparison.py](performance_comparison.py)

**Purpose**: Performance validation and benchmarking

**Key Features:**
- Golden record validation for numerical accuracy
- API overhead measurement
- SoA memory layout optimization
- OpenMP scaling analysis
- Integrator performance comparison
- Memory efficiency validation

**Run Command:**
```bash
python performance_comparison.py
```

**Prerequisites:**
- Multiple CPU cores (for OpenMP testing)
- Sufficient memory for large batches

**Expected Output:**
- Golden record validation results
- Performance benchmarks
- Scaling efficiency analysis
- Memory usage profiling

---

### 4. [scientific_analysis.py](scientific_analysis.py)

**Purpose**: Physics-focused analysis for fusion plasma research

**Key Features:**
- Orbit classification (trapped vs passing)
- Confinement analysis by flux surface
- Statistical analysis of particle behavior
- Parameter sweeps and convergence studies
- Physics visualization with matplotlib

**Run Command:**
```bash
python scientific_analysis.py
```

**Prerequisites:**
- matplotlib: `pip install matplotlib`
- scipy (optional): `pip install scipy`

**Expected Output:**
- Orbit classification statistics
- Confinement profiles
- Statistical analysis results
- Physics visualization plots

## Quick Start

1. **Build SIMPLE with Python support:**
   ```bash
   cd SIMPLE
   make  # Automatically includes Python API
   ```

2. **Download test data:**
   All examples automatically download the required VMEC equilibrium file.

3. **Run basic example:**
   ```bash
   cd python/examples
   python basic_batch_processing.py
   ```

## Example Workflow

For a complete scientific workflow:

```bash
# 1. Basic validation
python basic_batch_processing.py

# 2. Performance validation  
python performance_comparison.py

# 3. Large-scale simulation
python large_scale_streaming.py

# 4. Physics analysis
python scientific_analysis.py
```

## Common Issues and Solutions

### "pysimple module not available"

**Problem**: SIMPLE not built with Python support
**Solution**: 
```bash
cd SIMPLE
make clean
make  # Ensures Python support is included
```

### "Failed to download VMEC file"

**Problem**: Network connectivity or URL access
**Solution**: Download manually:
```bash
wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc
```

### Memory errors with large simulations

**Problem**: Insufficient memory for batch size
**Solution**: Use memory optimization:
```python
# Optimize batch size automatically
optimal_batch = simple.optimize_batch_size(
    target_memory_gb=4.0,  # Adjust for your system
    n_total=1_000_000
)
```

### Poor performance

**Problem**: Suboptimal thread configuration
**Solution**: Set OpenMP threads:
```python
results = simple.trace_orbits(
    particles,
    tmax=1000.0,
    openmp_threads=8  # Match your CPU cores
)
```

## Advanced Usage

### Custom Initialization

```python
# Initialize from existing data
particle_data = np.loadtxt("initial_conditions.dat")
particles = simple.ParticleBatch(len(particle_data))
particles.initialize_from_array(particle_data.T)  # Transpose to SoA
```

### Parameter Sweeps

```python
# Systematic parameter study
base_config = simple.create_configuration(
    vmec_file="wout.nc",
    tmax=1000.0,
    integrator='symplectic_midpoint'
)

dtau_values = [0.05, 0.1, 0.2, 0.5]
results = simple.parameter_sweep(
    base_config,
    'dtau',
    dtau_values,
    particles
)
```

### Streaming for Production

```python
# Production-scale simulation
stream_results = simple.process_large_simulation(
    vmec_file="wout.nc",
    n_total=10_000_000,  # 10 million particles
    tmax=2000.0,
    batch_size=100_000,
    output_file='production_results.h5',
    s_surface=0.9,
    memory_limit_gb=16.0
)
```

## Integration with Scientific Python

### NumPy Integration

```python
# Direct array operations
final_s = results.final_positions[0, :]
s_mean = np.mean(final_s)
s_std = np.std(final_s)

# Boolean indexing
confined_positions = results.final_positions[:, results.confined_mask]
```

### Matplotlib Visualization

```python
import matplotlib.pyplot as plt

# Poincaré plot
final_pos = results.final_positions
plt.scatter(final_pos[1, :], final_pos[2, :], s=1, alpha=0.5)
plt.xlabel('θ (poloidal)')
plt.ylabel('φ (toroidal)')
plt.show()
```

### HDF5 Analysis

```python
import h5py

# Efficient large dataset analysis
with h5py.File('results.h5', 'r') as f:
    for batch_name in f.keys():
        loss_times = f[batch_name]['loss_times'][:]
        confined_count = np.sum(loss_times == np.inf)
        print(f"{batch_name}: {confined_count} confined")
```

## Performance Guidelines

### Memory Optimization

- **Batch size**: Optimize for available memory using `optimize_batch_size()`
- **Streaming**: Use `process_large_simulation()` for >1M particles
- **Cleanup**: Delete large objects and call `gc.collect()` between batches

### CPU Optimization

- **OpenMP threads**: Set to number of physical cores
- **SoA access**: Use column access (`positions[:, particle_idx]`) over row access
- **Integrator choice**: `symplectic_midpoint` offers best performance/accuracy trade-off

### I/O Optimization

- **HDF5**: Use for large datasets with compression
- **NumPy**: Use for smaller datasets and interoperability
- **Batch processing**: Process in chunks rather than loading all data

## Next Steps

1. **Read the documentation**: `../docs/user_guide.md` and `../docs/api_reference.md`
2. **Explore your physics**: Modify examples for your specific research questions
3. **Scale up**: Use streaming capabilities for production simulations
4. **Contribute**: Report issues and suggest improvements

## Support

For questions and issues:

1. Check the troubleshooting section in `../docs/user_guide.md`
2. Review the API reference in `../docs/api_reference.md`
3. Examine the example source code for implementation details
4. File issues with specific error messages and system details

The SIMPLE Python API provides production-ready HPC capabilities while maintaining the scientific accuracy and performance of the underlying Fortran implementation.