"""
SIMPLE Python API - Pure Interface to Existing Fortran Functionality

Zero-copy interface layer providing direct access to proven SIMPLE Fortran
implementations without reimplementing computational logic.

Interface Design:
- All sampling delegates to existing samplers.f90 functions
- All simulation delegates to existing simple_main.f90 execution  
- All results access existing Fortran arrays via zero-copy wrappers
- <5% overhead vs direct Fortran execution

Core Principle: Python provides only parameter passing and result collection.
All computation performed by existing, proven Fortran implementations.
"""

from .core.batch import ParticleBatch, create_surface_batch, create_volume_batch, load_batch_from_file
from .core.results import BatchResults, ConfinementStats
from .core.simulation import (
    trace_orbits, create_configuration, benchmark_performance, 
    validate_golden_record, quick_simulation, parameter_sweep
)
from .samplers import SurfaceSampler, VolumeSampler, FileSampler
from .utils.memory import (
    process_large_simulation, estimate_memory_usage,
    optimize_batch_size, MemoryMonitor, run_million_particle_simulation
)

__version__ = "1.0.0"
__all__ = [
    # Core interface to existing Fortran arrays
    "ParticleBatch",
    "create_surface_batch", 
    "create_volume_batch",
    "load_batch_from_file",
    
    # Results wrappers around existing Fortran arrays
    "BatchResults", 
    "ConfinementStats",
    
    # Simulation interface to existing simple_main.f90
    "trace_orbits",
    "create_configuration",
    "benchmark_performance",
    "validate_golden_record", 
    "quick_simulation",
    "parameter_sweep",
    
    # Samplers interface to existing samplers.f90
    "SurfaceSampler",
    "VolumeSampler", 
    "FileSampler",
    
    # Memory utilities for large-scale coordination
    "process_large_simulation",
    "estimate_memory_usage",
    "optimize_batch_size",
    "MemoryMonitor",
    "run_million_particle_simulation"
]