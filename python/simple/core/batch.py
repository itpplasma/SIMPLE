"""
Pure interface wrapper around existing Fortran zstart arrays.

Zero-copy access to existing SoA layout without independent computation.
All particle initialization delegated to proven samplers.f90 functions.
"""

import numpy as np
from typing import Dict, Any, Optional
from dataclasses import dataclass
from ..backends.fortran import get_backend, FortranArrayWrapper


@dataclass
class Coordinates:
    """Zero-copy structured access to existing Fortran coordinate arrays"""
    s: np.ndarray       # Flux surface coordinate
    theta: np.ndarray   # Poloidal angle
    phi: np.ndarray     # Toroidal angle
    v_par: np.ndarray   # Parallel velocity
    mu: np.ndarray      # Magnetic moment
    
    @property
    def n_particles(self) -> int:
        """Number of particles from array size"""
        return self.s.shape[0] if self.s.ndim == 1 else self.s.shape[-1]


class ParticleBatch:
    """
    Zero-copy wrapper around existing zstart(5, n_particles) arrays.
    
    Pure interface to existing Fortran SoA layout - no independent computation.
    All functionality delegates to proven Fortran implementations.
    """
    
    def __init__(self, n_particles: int):
        """
        Create wrapper around existing Fortran allocation.
        
        Args:
            n_particles: Number of particles
        """
        if n_particles <= 0:
            raise ValueError("Number of particles must be positive")
        
        self.n_particles = n_particles
        self._backend = get_backend()
        self._arrays = self._backend.allocate_particle_arrays(n_particles)
    
    @classmethod
    def from_fortran_arrays(cls, fortran_arrays: np.ndarray) -> 'ParticleBatch':
        """
        Create batch from existing Fortran arrays.
        
        Args:
            fortran_arrays: Existing zstart array from Fortran
            
        Returns:
            ParticleBatch: Wrapper around existing arrays
        """
        if fortran_arrays.shape[0] != 5:
            raise ValueError("Fortran arrays must have shape (5, n_particles)")
        
        n_particles = fortran_arrays.shape[1]
        batch = cls(n_particles)
        
        # Copy data into Fortran-allocated arrays
        batch.positions[:] = fortran_arrays
        
        return batch
    
    @property
    def positions(self) -> np.ndarray:
        """
        Zero-copy view of existing zstart array.
        
        Returns:
            np.ndarray: Shape (5, n_particles) SoA layout from Fortran
        """
        return self._arrays.get_zstart_view()
    
    @property
    def coordinates(self) -> Coordinates:
        """
        Structured access to existing coordinate arrays.
        
        Returns:
            Coordinates: Zero-copy access to coordinate components
        """
        pos = self.positions
        return Coordinates(
            s=pos[0, :],        # Flux surface coordinate
            theta=pos[1, :],    # Poloidal angle
            phi=pos[2, :],      # Toroidal angle
            v_par=pos[3, :],    # Parallel velocity
            mu=pos[4, :]        # Magnetic moment
        )
    
    def initialize_from_samplers(self, vmec_file: str, method: str = 'surface', **kwargs):
        """
        Initialize using existing samplers.f90 functions.
        
        Args:
            vmec_file: Path to VMEC equilibrium file
            method: Sampling method ('surface', 'volume', 'file')
            **kwargs: Method-specific parameters
        """
        if method == 'surface':
            from ..samplers import SurfaceSampler
            sampler = SurfaceSampler(vmec_file)
            particle_data = sampler.sample_surface_fieldline(self.n_particles)
            
        elif method == 'volume':
            from ..samplers import VolumeSampler
            sampler = VolumeSampler(vmec_file)
            s_inner = kwargs.get('s_inner', 0.1)
            s_outer = kwargs.get('s_outer', 0.9)
            particle_data = sampler.sample_volume_single(self.n_particles, s_inner, s_outer)
            
        elif method == 'file':
            from ..samplers import FileSampler
            sampler = FileSampler(vmec_file)
            filename = kwargs.get('filename', 'start.dat')
            particle_data = sampler.load_from_file(filename)
            
        else:
            raise ValueError(f"Unknown method '{method}'. Use 'surface', 'volume', or 'file'")
        
        # Copy Fortran-generated data into our arrays
        self.positions[:] = particle_data
    
    def initialize_from_array(self, particle_data: np.ndarray):
        """
        Initialize from existing array data.
        
        Args:
            particle_data: Array with shape (5, n_particles) or (n_particles, 5)
        """
        if particle_data.shape == (5, self.n_particles):
            # SoA format - direct copy
            self.positions[:] = particle_data
        elif particle_data.shape == (self.n_particles, 5):
            # AoS format - transpose to SoA
            self.positions[:] = particle_data.T
        else:
            raise ValueError(f"Invalid particle data shape {particle_data.shape}, "
                           f"expected (5, {self.n_particles}) or ({self.n_particles}, 5)")
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Export data for serialization/analysis.
        
        Returns:
            Dict: Serializable representation
        """
        return {
            'coordinates': self.positions.copy(),
            'n_particles': self.n_particles,
            'metadata': self._arrays.get_metadata()
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ParticleBatch':
        """
        Reconstruct from serialized data.
        
        Args:
            data: Dictionary from to_dict()
            
        Returns:
            ParticleBatch: Reconstructed batch
        """
        batch = cls(data['n_particles'])
        batch.positions[:] = data['coordinates']
        return batch
    
    def copy(self) -> 'ParticleBatch':
        """
        Create independent copy.
        
        Returns:
            ParticleBatch: Independent copy with own memory
        """
        new_batch = ParticleBatch(self.n_particles)
        new_batch.positions[:] = self.positions.copy()
        return new_batch
    
    def slice(self, start: int, end: Optional[int] = None) -> 'ParticleBatch':
        """
        Create subset batch.
        
        Args:
            start: Starting particle index
            end: Ending particle index (exclusive)
            
        Returns:
            ParticleBatch: New batch with subset
        """
        if end is None:
            end = self.n_particles
        
        if not (0 <= start < end <= self.n_particles):
            raise ValueError(f"Invalid slice range [{start}:{end}] for {self.n_particles} particles")
        
        subset_size = end - start
        subset_batch = ParticleBatch(subset_size)
        subset_batch.positions[:] = self.positions[:, start:end].copy()
        
        return subset_batch
    
    def get_particle(self, index: int) -> np.ndarray:
        """
        Get single particle coordinates (cache-friendly column access).
        
        Args:
            index: Particle index
            
        Returns:
            np.ndarray: Shape (5,) with [s, theta, phi, v_par, mu]
        """
        if not (0 <= index < self.n_particles):
            raise IndexError(f"Particle index {index} out of range [0, {self.n_particles})")
        
        return self.positions[:, index].copy()
    
    def set_particle(self, index: int, coordinates: np.ndarray):
        """
        Set single particle coordinates.
        
        Args:
            index: Particle index
            coordinates: Array with shape (5,)
        """
        if not (0 <= index < self.n_particles):
            raise IndexError(f"Particle index {index} out of range [0, {self.n_particles})")
        
        if coordinates.shape != (5,):
            raise ValueError(f"Coordinates must have shape (5,), got {coordinates.shape}")
        
        self.positions[:, index] = coordinates
    
    def validate_layout(self) -> Dict[str, Any]:
        """
        Validate SoA memory layout.
        
        Returns:
            Dict: Layout validation results
        """
        positions = self.positions
        
        return {
            'shape': positions.shape,
            'c_contiguous': positions.flags.c_contiguous,
            'strides': getattr(positions, 'strides', None),
            'dtype': positions.dtype,
            'fortran_allocated': True  # Always true - from backend
        }
    
    def __len__(self) -> int:
        """Return number of particles"""
        return self.n_particles
    
    def __repr__(self) -> str:
        """String representation"""
        return f"ParticleBatch(n_particles={self.n_particles})"
    
    def __str__(self) -> str:
        """Detailed string representation"""
        coords = self.coordinates
        s_range = f"[{coords.s.min():.3f}, {coords.s.max():.3f}]"
        return (f"ParticleBatch({self.n_particles} particles)\n"
                f"  s range: {s_range}\n"
                f"  SoA shape: {self.positions.shape}")


# Convenience functions for common initialization patterns
def create_surface_batch(n_particles: int, vmec_file: str, s: float = 0.9) -> ParticleBatch:
    """
    Create batch using existing surface sampling.
    
    Args:
        n_particles: Number of particles
        vmec_file: VMEC equilibrium file
        s: Flux surface coordinate
        
    Returns:
        ParticleBatch: Initialized batch
    """
    batch = ParticleBatch(n_particles)
    batch.initialize_from_samplers(vmec_file, method='surface', s=s)
    return batch


def create_volume_batch(
    n_particles: int, 
    vmec_file: str, 
    s_inner: float = 0.1, 
    s_outer: float = 0.9
) -> ParticleBatch:
    """
    Create batch using existing volume sampling.
    
    Args:
        n_particles: Number of particles
        vmec_file: VMEC equilibrium file
        s_inner: Inner flux surface
        s_outer: Outer flux surface
        
    Returns:
        ParticleBatch: Initialized batch
    """
    batch = ParticleBatch(n_particles)
    batch.initialize_from_samplers(vmec_file, method='volume', s_inner=s_inner, s_outer=s_outer)
    return batch


def load_batch_from_file(vmec_file: str, filename: str) -> ParticleBatch:
    """
    Load batch using existing file loading.
    
    Args:
        vmec_file: VMEC equilibrium file
        filename: File containing particle data
        
    Returns:
        ParticleBatch: Loaded batch
    """
    from ..samplers import FileSampler
    sampler = FileSampler(vmec_file)
    particle_data = sampler.load_from_file(filename)
    
    return ParticleBatch.from_fortran_arrays(particle_data)