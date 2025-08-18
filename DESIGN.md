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

### 1. **Primary Use Case Optimization**
```python
# OPTIMIZED FOR: Massive parallel orbit tracing
config = {
    'vmec_file': 'wout.nc',
    'nparticles': 10000,      # Thousands of particles
    'tmax': 1000.0,           # Long simulation times
    'field_type': 'vmec_boozer'
}

sim = simple.Simulation(config)
results = sim.run()  # OpenMP parallelized, minimal Python overhead
```

### 2. **Minimal Overhead Requirement**
- Python API must add <5% overhead to Fortran execution
- Configuration conversion cached for repeated use
- Bulk data transfer optimization
- Memory management aligned with Fortran patterns

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

### 3. **Performance Validation**
- Benchmark testing at each phase
- Memory usage profiling
- OpenMP scaling verification
- Scientific accuracy validation

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