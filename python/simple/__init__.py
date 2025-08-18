"""
SIMPLE Python API - Batch-Oriented HPC Interface

High-performance Python interface for SIMPLE particle orbit tracing
with zero-copy access to existing Fortran SoA data structures.

API Design:
- ParticleBatch: Zero-copy wrapper around existing zstart(5, n_particles) arrays
- BatchResults: Performance-optimized access to result arrays
- trace_orbits(): Batch processing leveraging existing OpenMP parallelization

Performance Guarantees:
- <5% overhead vs direct Fortran execution
- Zero-copy access to existing SoA structures  
- Memory-efficient streaming for millions of particles
- GPU-ready architecture with existing SoA layout
"""

from .core.batch import ParticleBatch
from .core.results import BatchResults, ConfinementStats
from .core.simulation import (
    trace_orbits, create_configuration, benchmark_performance, 
    validate_golden_record, quick_simulation, parameter_sweep
)
from .samplers import SurfaceSampler, VolumeSampler, FileSampler
from .utils.memory import (
    ParticleBatchStream, process_large_simulation, estimate_memory_usage,
    optimize_batch_size, MemoryMonitor
)

__version__ = "1.0.0"
__all__ = [
    "ParticleBatch",
    "BatchResults", 
    "ConfinementStats",
    "trace_orbits",
    "create_configuration",
    "benchmark_performance",
    "validate_golden_record", 
    "quick_simulation",
    "parameter_sweep",
    "SurfaceSampler",
    "VolumeSampler", 
    "FileSampler",
    "ParticleBatchStream",
    "process_large_simulation",
    "estimate_memory_usage",
    "optimize_batch_size",
    "MemoryMonitor"
]