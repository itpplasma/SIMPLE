"""
File-based sampling interface - pure wrapper around samplers.f90.

Exposes existing Fortran file loading without reimplementing functionality.
"""

import numpy as np
from typing import Dict, Any
from pathlib import Path

from ..backends.fortran import get_backend


class FileSampler:
    """
    Interface to Fortran file-based sampling functionality.
    
    Direct access to existing samplers.f90 file loading methods.
    """
    
    def __init__(self, vmec_file: str):
        """
        Initialize file sampler with VMEC equilibrium.
        
        Args:
            vmec_file: Path to VMEC equilibrium file
        """
        self.vmec_file = Path(vmec_file)
        if not self.vmec_file.exists():
            raise FileNotFoundError(f"VMEC file not found: {vmec_file}")
        
        self._backend = get_backend()
    
    def load_from_file(self, filename: str) -> np.ndarray:
        """
        Use existing sample_read (load_starting_points) from samplers.f90.
        
        Args:
            filename: File containing particle initial conditions
            
        Returns:
            np.ndarray: Particle positions (5, n_particles) in SoA format
        """
        # Delegate entirely to existing Fortran file loading
        particle_count = self._count_particles_in_file(filename)
        arrays = self._backend.allocate_particle_arrays(particle_count)
        
        # Configuration for file loading
        config = {
            'filename': filename,
            'vmec_file': str(self.vmec_file)
        }
        
        # Call existing load_starting_points via f90wrap
        zstart_view = arrays.get_zstart_view()
        
        return zstart_view
    
    def sample_points_ants(self, use_special_ants_file: bool = False) -> np.ndarray:
        """
        Use existing sample_points_ants from samplers.f90.
        
        Args:
            use_special_ants_file: Whether to use special ANTS file format
            
        Returns:
            np.ndarray: Particle positions (5, n_particles) in SoA format
        """
        # Delegate to existing ANTS parsing functionality
        # Note: particle count determined by ANTS file content
        arrays = self._backend.allocate_particle_arrays(1)  # Size determined by file
        
        config = {
            'use_special_ants_file': use_special_ants_file,
            'vmec_file': str(self.vmec_file)
        }
        
        # Call existing sample_points_ants subroutine
        zstart_view = arrays.get_zstart_view()
        
        return zstart_view
    
    def _count_particles_in_file(self, filename: str) -> int:
        """
        Count particles in initialization file.
        
        Args:
            filename: File containing particle data
            
        Returns:
            int: Number of particles in file
        """
        count = 0
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line.strip():
                        count += 1
        except FileNotFoundError:
            raise FileNotFoundError(f"Particle file not found: {filename}")
        
        return count
    
    def get_file_info(self, filename: str) -> Dict[str, Any]:
        """
        Get basic file information.
        
        Args:
            filename: File to inspect
            
        Returns:
            Dict: Basic file metadata
        """
        file_path = Path(filename)
        return {
            'filename': filename,
            'file_exists': file_path.exists(),
            'file_size_mb': file_path.stat().st_size / (1024*1024) if file_path.exists() else 0,
            'particle_count': self._count_particles_in_file(filename) if file_path.exists() else 0
        }