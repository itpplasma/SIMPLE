"""
Backend implementations for SIMPLE Python API.

Modules:
- fortran: f90wrap integration with existing pysimple module
- native: Future direct memory interface implementation
"""

from .fortran import FortranBackend, FortranArrayWrapper, FortranResultWrapper

__all__ = [
    "FortranBackend",
    "FortranArrayWrapper", 
    "FortranResultWrapper"
]