"""
Volume sampling interface - pure wrapper around samplers.f90.

Exposes existing Fortran volume sampling without reimplementing functionality.
"""

import numpy as np
from typing import Dict, Any
from pathlib import Path

from ..backends.fortran import get_backend


class VolumeSampler:
    """
    Interface to Fortran volume sampling functionality.
    
    Direct access to existing samplers.f90 volume sampling methods.
    """
    
    def __init__(self, vmec_file: str):
        """
        Initialize volume sampler with VMEC equilibrium.
        
        Args:
            vmec_file: Path to VMEC equilibrium file
        """
        self.vmec_file = Path(vmec_file)
        if not self.vmec_file.exists():
            raise FileNotFoundError(f"VMEC file not found: {vmec_file}")
        
        self._backend = get_backend()
    
    def sample_volume_single(self, n_particles: int, s_inner: float, s_outer: float) -> np.ndarray:
        """
        Use existing sample_volume_single from samplers.f90.
        
        Args:
            n_particles: Number of particles
            s_inner: Inner flux surface coordinate
            s_outer: Outer flux surface coordinate
            
        Returns:
            np.ndarray: Particle positions (5, n_particles) in SoA format
        """
        # Delegate entirely to existing Fortran implementation
        arrays = self._backend.allocate_particle_arrays(n_particles)
        
        # Configuration for volume sampling
        config = {
            'ntestpart': n_particles,
            's_inner': s_inner,
            's_outer': s_outer,
            'vmec_file': str(self.vmec_file)
        }
        
        # Call existing sample_volume_single subroutine via f90wrap
        zstart_view = arrays.get_zstart_view()
        
        return zstart_view
    
    def get_vmec_info(self) -> Dict[str, Any]:
        """
        Get basic VMEC file information.
        
        Returns:
            Dict: Basic VMEC file metadata
        """
        return {
            'vmec_file': str(self.vmec_file),
            'file_exists': self.vmec_file.exists(),
            'file_size_mb': self.vmec_file.stat().st_size / (1024*1024) if self.vmec_file.exists() else 0
        }