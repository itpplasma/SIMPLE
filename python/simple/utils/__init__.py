"""
Utility modules for SIMPLE Python API.

Modules:
- memory: Memory-efficient streaming and large dataset utilities
- visualization: Batch visualization tools (future)
"""

from .memory import (
    process_large_simulation, estimate_memory_usage,
    optimize_batch_size, MemoryMonitor, run_million_particle_simulation
)

__all__ = [
    "process_large_simulation",
    "estimate_memory_usage", 
    "optimize_batch_size",
    "MemoryMonitor",
    "run_million_particle_simulation"
]