# SIMPLE Design Document

## Strategic Vision: Python-First API Design

**CORE PRINCIPLE**: Python serves as **performance-first API prototype** that demonstrates clean design patterns for subsequent Fortran API implementation, with eventual API convergence.

### Four-Phase Strategic Roadmap

#### **Phase 1**: Python API Prototype (`simple` Module)
- **PRIMARY USE CASE**: OpenMP parallelized massive orbit tracer (performance-first)
- **MODULE NAME**: `simple` (not `pysimple` or `simple_api`)
- **NAMING CONVENTION**: Python uses PEP8 conventions (`Config`, `Simulation`, `Orbit`)
- **PURPOSE**: Demonstrate clean API patterns that will be implemented in Fortran

#### **Phase 2**: Python-Fortran Name Translation Layer  
- **AUTOMATIC TRANSLATION**: `typename_t` ↔ `PEP8` convention mapping
- **TRANSPARENT**: Users see only their preferred convention
- **SEAMLESS INTEGRATION**: Clean interface between Python and Fortran

#### **Phase 3**: Fortran API Refactoring
- **PATTERN IMPLEMENTATION**: Apply proven Python patterns to Fortran
- **NAMING CONVENTION**: Fortran uses `typename_t` convention (`config_t`, `simulation_t`, `orbit_t`)
- **ATOMIC STEPS**: Safe incremental refactoring preserving `simple.x` functionality
- **GOLDEN RECORD**: Validation ensures no regressions

#### **Phase 4**: Python as Thin Wrapper
- **END STATE**: Python becomes lightweight wrapper over clean Fortran API
- **PERFORMANCE**: Remains in Fortran, Python provides convenience
- **API CONVERGENCE**: Both interfaces expose same clean patterns

## Current Architecture Analysis

### Strengths
- **Performance**: Highly optimized symplectic integrators with OpenMP parallelization
- **Scientific Accuracy**: Proven algorithms for stellarator particle tracing
- **Robustness**: Stable codebase with extensive real-world usage

### API Modernization Needs

#### 1. **Global State Dependencies**
```fortran
! Current: Global module variables
module params
    real(dp) :: tmax, dtau
    integer :: nparticles, field_type
end module

! Target: Encapsulated configuration
type :: config_t
    real(dp) :: tmax, dtau
    integer :: nparticles, field_type
contains
    procedure :: from_namelist
    procedure :: validate
end type
```

#### 2. **Rigid Initialization Order**
```fortran
! Current: Rigid sequence with side effects
call read_config()
call init_field()  
call params_init()
call run()

! Target: Clean encapsulated interface
type(simulation_t) :: sim
sim = create_simulation(config)
results = sim%run()
```

#### 3. **Low-Level API Exposure**
```python
# Current: Complex multi-step setup
import pysimple
pysimple.read_config('simple.in')
pysimple.init_field()
pysimple.params_init()
pysimple.run()

# Target: Clean high-level interface  
import simple
config = {'vmec_file': 'wout.nc', 'nparticles': 10000}
sim = simple.Simulation(config)
results = sim.run()
```

## Design Patterns and Conventions

### Naming Translation Strategy

#### Python Convention (PEP8)
```python
class Config:
    def from_dict(self, config_dict): pass
    def to_namelist_string(self): pass

class Simulation:
    def __init__(self, config): pass
    def run(self): pass

class Orbit:
    def trace(self, tmax): pass
```

#### Fortran Convention (`typename_t`)
```fortran
type :: config_t
contains
    procedure :: from_dict => config_from_dict
    procedure :: to_namelist_string => config_to_namelist
end type

type :: simulation_t  
contains
    procedure :: init => simulation_init
    procedure :: run => simulation_run
end type

type :: orbit_t
contains
    procedure :: trace => orbit_trace
end type
```

#### Automatic Translation Layer
```python
# Name translator handles convention mapping
class NameTranslator:
    def python_to_fortran(self, name):
        return f"{name.lower()}_t"
    
    def method_python_to_fortran(self, class_name, method_name):
        return f"{class_name.lower()}_{method_name}"
```

### Core Module Architecture

#### Python Module Structure (`simple/`)
```
simple/
├── __init__.py          # Public API (Config, Simulation, Orbit)
├── config.py           # Configuration management prototype
├── simulation.py       # Batch simulation interface
├── orbit.py            # Single orbit analysis
├── translation.py      # Name convention translator
├── _backend.py         # f90wrap interface to Fortran
└── visualization.py    # Built-in plotting
```

#### Fortran Module Structure (Target)
```
src/api/
├── simple_config.f90      # config_t derived type
├── simple_simulation.f90  # simulation_t derived type  
├── simple_orbit.f90       # orbit_t derived type
├── simple_api.f90         # High-level API functions
└── simple_translation.f90 # Name translation utilities
```

## Performance-First Design Principles

**CRITICAL INSIGHT**: Current Fortran implementation is ALREADY HPC/SoA optimized:
- `zstart(1:5, ntestpart)` - Structure of Arrays pattern ✓
- OpenMP particle parallelization ✓  
- Cache-friendly column access ✓
- Million+ particle capability ✓

**DESIGN GOAL**: Create batch-oriented APIs around existing high-performance core.

### 1. **Batch-Oriented API Design**
```python
# LEVERAGES EXISTING: SoA data structures, OpenMP parallelization
import simple

# Batch processing (primary use case)
particles = simple.ParticleBatch.from_vmec_surface(
    vmec_file='wout.nc',
    n_particles=100_000,  # Leverages existing SoA: zstart(5, n_particles)
    surface_s=0.9
)

# GPU-ready batch operations
results = simple.trace_orbits(
    particles,
    tmax=1000.0,
    integrator='symplectic_midpoint'  # Uses existing optimized integrators
)

# Direct access to underlying SoA arrays
positions = results.final_positions    # Shape: (5, n_particles) - existing zstart format
loss_times = results.loss_times        # Shape: (n_particles,) - existing times_lost format
```

### 2. **GPU-Ready SoA Architecture Strategy**
```python
# LEVERAGE EXISTING: SoA layout perfect for GPU acceleration
class ParticleBatch:
    def to_gpu(self) -> 'GPUParticleBatch':
        # Existing SoA layout transfers efficiently to GPU
        gpu_positions = cupy.asarray(self.positions)  # Shape: (5, n_particles)
        return GPUParticleBatch(gpu_positions)
    
    def from_gpu(self) -> 'ParticleBatch':
        # Efficient GPU->CPU transfer maintaining SoA layout
        cpu_positions = cupy.asnumpy(self._gpu_positions)
        return ParticleBatch(cpu_positions)

# Future GPU kernel interface (leverages existing SoA)
def trace_orbits_gpu(particles: GPUParticleBatch, **kwargs):
    # GPU kernels operate on existing SoA memory layout
    # positions[5, n_particles] - coalesced memory access
    # One thread per particle - maps to existing OpenMP pattern
    pass
```

### 3. **Batch Processing Optimization** 
```python
# STREAM PROCESSING: Handle memory constraints efficiently
def process_large_simulation(n_total: int, batch_size: int = 100_000):
    """Process millions of particles in memory-efficient batches"""
    results = []
    
    for start_idx in range(0, n_total, batch_size):
        end_idx = min(start_idx + batch_size, n_total)
        batch_size_actual = end_idx - start_idx
        
        # Create batch using existing SoA allocation patterns
        batch = simple.ParticleBatch(n_particles=batch_size_actual)
        batch.initialize_surface(vmec_file='wout.nc', s=0.9)
        
        # Use existing OpenMP parallelization
        batch_results = simple.trace_orbits(batch, tmax=1000.0)
        results.append(batch_results)
        
        # Clean memory between batches
        del batch, batch_results
    
    return combine_batch_results(results)
```

### 4. **Minimal Overhead Requirement**
- Python API must add <5% overhead to Fortran execution
- Zero-copy array access to existing SoA structures
- Configuration conversion cached for repeated use
- Memory management aligned with existing Fortran patterns

### 3. **Scientific Computing Integration**
```python
# Clean integration with scientific Python ecosystem
import simple
import numpy as np
import matplotlib.pyplot as plt

config = simple.Config.from_file('simple.in')
config.nparticles = 5000

sim = simple.Simulation(config)
results = sim.run()

# Direct NumPy array access to results
confined_fraction = results.confined_fraction  # NumPy array
loss_times = results.loss_times              # NumPy array

# Built-in visualization
results.plot_confinement()
results.plot_poincare_section()
```

## Implementation Quality Standards

### 1. **Atomic Refactoring Principle**
Each implementation step must:
- Preserve existing `simple.x` functionality exactly
- Pass all golden record tests unchanged
- Allow rollback if issues discovered
- Be reviewable in isolation

### 2. **Convention Consistency**
- **Fortran**: Strict `typename_t` convention throughout
- **Python**: Strict PEP8 convention throughout  
- **Translation**: Automatic and transparent
- **Documentation**: Consistent with target convention

### 3. **Performance Validation: HPC Characteristics Preservation**

#### **SoA Layout Validation**
```python
# Verify existing SoA patterns are preserved and accessible
def validate_soa_performance():
    import simple
    import numpy as np
    import time
    
    # Test 1: SoA memory layout efficiency
    n_particles = 1_000_000
    batch = simple.ParticleBatch(n_particles)
    
    # Verify contiguous memory layout (existing pattern)
    positions = batch.positions  # Shape: (5, n_particles)
    assert positions.flags.c_contiguous, "SoA layout must be contiguous"
    
    # Test 2: Column access performance (existing cache pattern)
    start_time = time.time()
    for i in range(min(1000, n_particles)):
        particle_coords = positions[:, i]  # Cache-friendly column access
    column_access_time = time.time() - start_time
    
    # Test 3: Row access performance (should be slower)
    start_time = time.time()
    for j in range(5):
        coord_array = positions[j, :1000]  # Row access
    row_access_time = time.time() - start_time
    
    # Validate existing cache-friendly pattern is preserved
    assert column_access_time < row_access_time * 0.5, "Column access should be faster (SoA benefit)"

# Test 4: OpenMP scaling validation
def validate_openmp_scaling():
    """Verify existing OpenMP performance is preserved"""
    import simple
    import os
    
    results_1_thread = {}
    results_8_threads = {}
    
    # Test with 1 thread
    os.environ['OMP_NUM_THREADS'] = '1'
    start_time = time.time()
    result_1 = simple.trace_orbits(particles, tmax=100.0)
    results_1_thread['time'] = time.time() - start_time
    
    # Test with 8 threads 
    os.environ['OMP_NUM_THREADS'] = '8'
    start_time = time.time()
    result_8 = simple.trace_orbits(particles, tmax=100.0)
    results_8_threads['time'] = time.time() - start_time
    
    # Verify scaling efficiency (existing OpenMP patterns)
    scaling_efficiency = results_1_thread['time'] / (8 * results_8_threads['time'])
    assert scaling_efficiency > 0.6, f"OpenMP scaling too low: {scaling_efficiency}"
    
    # Verify results are identical (deterministic with existing algorithms)
    np.testing.assert_array_equal(result_1.loss_times, result_8.loss_times)

# Test 5: Memory efficiency validation
def validate_memory_efficiency():
    """Ensure no memory overhead vs existing implementation"""
    import psutil
    import simple
    
    # Baseline memory usage
    process = psutil.Process()
    baseline_memory = process.memory_info().rss
    
    # Create batch using new API
    n_particles = 100_000
    batch = simple.ParticleBatch(n_particles)
    api_memory = process.memory_info().rss
    
    # Expected memory for SoA arrays: 5 * n_particles * 8 bytes + additional arrays
    expected_particle_memory = 5 * n_particles * 8  # positions
    expected_particle_memory += n_particles * 8     # loss_times
    expected_particle_memory += n_particles * 8     # trap_par
    expected_particle_memory += n_particles * 8     # perp_inv
    
    actual_memory_increase = api_memory - baseline_memory
    overhead_ratio = actual_memory_increase / expected_particle_memory
    
    assert overhead_ratio < 1.1, f"Memory overhead too high: {overhead_ratio:.2f}"
```

#### **Golden Record Performance Testing**
```python
def validate_golden_record_performance():
    """Compare against existing Fortran implementation performance"""
    
    # Performance benchmarks must match existing simple.x within 5%
    fortran_time = benchmark_fortran_simple()  # Existing implementation
    python_time = benchmark_python_api()       # New API
    
    performance_ratio = python_time / fortran_time
    assert performance_ratio < 1.05, f"Python API too slow: {performance_ratio:.3f}x"
    
    # Results must be bit-identical
    fortran_results = load_fortran_results()
    python_results = run_python_equivalent() 
    
    np.testing.assert_array_equal(
        fortran_results['times_lost'], 
        python_results.loss_times,
        "Results must be identical to existing implementation"
    )

def benchmark_memory_scaling():
    """Validate memory usage scales linearly with particle count"""
    particle_counts = [10_000, 100_000, 1_000_000]
    memory_usage = []
    
    for n in particle_counts:
        process = psutil.Process()
        before = process.memory_info().rss
        
        batch = simple.ParticleBatch(n)
        batch.initialize_surface('wout.nc', s=0.9)
        
        after = process.memory_info().rss
        memory_usage.append((after - before) / n)  # Per-particle memory
        
        del batch
    
    # Memory per particle should be constant (linear scaling)
    assert all(abs(m - memory_usage[0]) / memory_usage[0] < 0.1 for m in memory_usage), \
        "Memory scaling must be linear"
```

#### **GPU Preparation Validation**
```python
def validate_gpu_readiness():
    """Ensure SoA layout is optimal for future GPU acceleration"""
    import simple
    
    batch = simple.ParticleBatch(100_000)
    positions = batch.positions  # Shape: (5, n_particles)
    
    # Test 1: Contiguous memory for efficient GPU transfer
    assert positions.flags.c_contiguous, "Must be C-contiguous for GPU transfer"
    
    # Test 2: Coalesced memory access pattern
    # GPU threads would access positions[:, thread_id] - coalesced
    for thread_block in range(0, 1000, 256):  # Simulate GPU thread blocks
        thread_ids = range(thread_block, min(thread_block + 256, 1000))
        particle_data = positions[:, thread_ids]  # Coalesced access pattern
        assert particle_data.shape[0] == 5, "Coordinate access must be contiguous"
    
    # Test 3: Memory layout matches GPU optimization requirements
    stride_coords = positions.strides[0]  # Stride between coordinates
    stride_particles = positions.strides[1]  # Stride between particles  
    
    # For coalesced access, particle stride should be smaller
    assert stride_particles == 8, "Particle stride optimal for GPU coalescing" 
    assert stride_coords == 8 * positions.shape[1], "Coordinate stride correct for SoA"
```

## Success Metrics

### Phase 1 (Python Prototype)
- [ ] Performance: <5% overhead vs direct Fortran
- [ ] Usability: Single orbit in <5 lines of code
- [ ] Isolation: Multiple simulations without interference
- [ ] Golden Record: Exact equivalence with existing implementation

### Phase 2 (Name Translation)
- [ ] Transparent: Users unaware of translation layer
- [ ] Complete: All API elements correctly mapped
- [ ] Consistent: Predictable naming patterns
- [ ] Performance: No measurable translation overhead

### Phase 3 (Fortran Refactoring)  
- [ ] API Parity: Fortran matches Python interface patterns
- [ ] Convention: Consistent `typename_t` throughout
- [ ] Compatibility: `simple.x` unchanged externally
- [ ] Quality: Modern Fortran patterns applied

### Phase 4 (Convergence)
- [ ] Thin Wrapper: Python primarily delegates to Fortran
- [ ] Performance: Near-native Fortran performance through Python
- [ ] Maintenance: Single source of truth for API patterns
- [ ] User Experience: Consistent interface across languages

## Benefits of This Approach

### 1. **Risk Mitigation**
- Python prototype validates patterns before Fortran changes
- Atomic refactoring allows safe rollback
- Golden record testing prevents regressions
- User feedback guides API design

### 2. **Development Efficiency**
- Python enables rapid API iteration
- Patterns proven before costly Fortran implementation
- Incremental value delivery throughout process
- Clear blueprint for Fortran modernization

### 3. **Long-term Maintainability**  
- Modern API patterns in both languages
- Consistent naming conventions
- Clean separation of concerns
- Foundation for future enhancements

This design creates a pathway from the current Fortran-centric codebase to a modern, dual-language API while preserving performance and maintaining backward compatibility throughout the transition.