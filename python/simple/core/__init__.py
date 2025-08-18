"""
Core API components for batch-oriented particle processing.

Modules:
- batch: ParticleBatch class with zero-copy SoA access
- results: BatchResults class with performance methods
- simulation: High-level simulation interface and configuration
"""

from .batch import ParticleBatch, Coordinates
from .results import BatchResults, ConfinementStats
from .simulation import trace_orbits, create_configuration

__all__ = [
    "ParticleBatch",
    "Coordinates", 
    "BatchResults",
    "ConfinementStats",
    "trace_orbits",
    "create_configuration"
]