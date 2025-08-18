# SIMPLE Python API Examples

This directory contains minimal examples demonstrating the pure interface layer to existing Fortran functionality.

## Design Philosophy

The Python API is a **pure interface layer** that provides zero-copy access to existing SIMPLE Fortran implementations without reimplementing any computational logic.

All functionality delegates to proven Fortran code:
- Particle sampling → `samplers.f90`
- Orbit integration → `simple_main.f90` 
- Field evaluation → `field_*.f90`
- Results processing → existing Fortran arrays

## Basic Usage

### Quick Start

```python
import simple

# 1. Initialize particles using existing samplers.f90
sampler = simple.SurfaceSampler("wout.nc")
particle_data = sampler.sample_surface_fieldline(n_particles=1000)
particles = simple.ParticleBatch.from_fortran_arrays(particle_data)

# 2. Run simulation using existing simple_main.f90
results = simple.trace_orbits(
    particles,
    tmax=500.0,
    integrator='symplectic_midpoint'
)

# 3. Access results through zero-copy Fortran array wrappers
stats = results.confinement_statistics()
print(f"Confined fraction: {stats.confined_fraction:.3f}")
```

## Examples

### basic_usage.py
Comprehensive example showing:
- Surface, volume, and file-based particle initialization
- Different integrator methods
- Results access patterns
- Pure interface usage patterns

**Run**: `python basic_usage.py`

## Requirements

1. **Build SIMPLE with Python support**:
   ```bash
   cd /path/to/SIMPLE
   make  # Builds both Fortran and Python interfaces
   ```

2. **Test VMEC file**:
   Download automatically or manually:
   ```bash
   wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc
   ```

## Testing

Comprehensive testing is in `test/python/`:
- `test_batch_api.py` - API functionality tests
- `test_performance_validation.py` - Performance overhead validation  
- `test_vmec.py` - VMEC integration tests

Run tests:
```bash
cd test/python
python -m pytest -v
```

## Zero-Copy Interface Design

The Python API exposes existing Fortran functionality through:

1. **Samplers**: Direct access to `samplers.f90` functions
   - `SurfaceSampler` → `sample_surface_fieldline()`
   - `VolumeSampler` → `sample_volume_single()`
   - `FileSampler` → `sample_read()`, `sample_points_ants()`

2. **Execution**: Delegate to `simple_main.f90`
   - `trace_orbits()` → calls existing `run()` subroutine
   - Zero computational logic in Python

3. **Results**: Zero-copy wrappers around Fortran arrays
   - `BatchResults` → wraps `times_lost`, `zend`, etc.
   - Direct access to existing SoA layout

## Performance

The interface provides <5% overhead vs direct Fortran execution by:
- Zero-copy access to existing `zstart(5, n_particles)` arrays
- Direct delegation to existing OpenMP parallelization
- No independent computation in Python layer

## API Reference

See `python/docs/api_reference.md` for detailed documentation.