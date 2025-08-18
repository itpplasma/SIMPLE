"""
Surface sampling for particle initialization on flux surfaces.

Integrates with existing samplers.f90 functionality for consistent 
particle initialization patterns across Python and Fortran interfaces.
"""

import numpy as np
from typing import Optional, Dict, Any
from pathlib import Path

from ..core.batch import ParticleBatch


class SurfaceSampler:
    """
    Surface sampling for particle initialization on flux surfaces.
    
    This sampler provides various methods for distributing particles
    on magnetic flux surfaces using existing SIMPLE sampling algorithms.
    """
    
    def __init__(self, vmec_file: str):
        """
        Initialize surface sampler with VMEC equilibrium.
        
        Args:
            vmec_file: Path to VMEC equilibrium file
        """
        self.vmec_file = Path(vmec_file)
        if not self.vmec_file.exists():
            raise FileNotFoundError(f"VMEC file not found: {vmec_file}")
        
        self._cached_equilibrium = None
    
    def uniform_poloidal(
        self,
        n_particles: int,
        s: float,
        velocity_spread: float = 0.1,
        pitch_angle_range: tuple = (-1.0, 1.0),
        seed: Optional[int] = None
    ) -> ParticleBatch:
        """
        Uniform distribution in poloidal angle on flux surface.
        
        Args:
            n_particles: Number of particles to generate
            s: Flux surface coordinate (0 < s < 1)
            velocity_spread: Relative velocity spread around unity
            pitch_angle_range: Range for pitch angle (v_parallel/v)
            seed: Random seed for reproducibility
            
        Returns:
            ParticleBatch: Initialized particle batch
        """
        if not (0 < s < 1):
            raise ValueError("Flux surface s must be between 0 and 1")
        
        if seed is not None:
            np.random.seed(seed)
        
        batch = ParticleBatch(n_particles)
        positions = batch.positions
        
        # Flux surface coordinate with small random variation
        positions[0, :] = s + velocity_spread * 0.1 * (np.random.random(n_particles) - 0.5)
        
        # Uniform poloidal distribution
        positions[1, :] = np.random.uniform(0, 2*np.pi, n_particles)
        
        # Random toroidal distribution
        positions[2, :] = np.random.uniform(0, 2*np.pi, n_particles)
        
        # Velocity with spread
        base_velocity = 1.0
        velocity_variation = 1.0 + velocity_spread * (np.random.random(n_particles) - 0.5)
        positions[3, :] = base_velocity * velocity_variation
        
        # Pitch angle (v_parallel / v)
        pitch_min, pitch_max = pitch_angle_range
        positions[4, :] = np.random.uniform(pitch_min, pitch_max, n_particles)
        
        return batch
    
    def flux_surface_grid(
        self,
        n_theta: int,
        n_phi: int,
        s: float,
        velocity: float = 1.0,
        pitch_angles: Optional[np.ndarray] = None
    ) -> ParticleBatch:
        """
        Regular grid on flux surface in (theta, phi) coordinates.
        
        Args:
            n_theta: Number of poloidal grid points
            n_phi: Number of toroidal grid points  
            s: Flux surface coordinate
            velocity: Particle velocity magnitude
            pitch_angles: Array of pitch angles (or uniform if None)
            
        Returns:
            ParticleBatch: Grid-initialized particle batch
        """
        if not (0 < s < 1):
            raise ValueError("Flux surface s must be between 0 and 1")
        
        n_particles = n_theta * n_phi
        batch = ParticleBatch(n_particles)
        positions = batch.positions
        
        # Create regular grid
        theta_grid = np.linspace(0, 2*np.pi, n_theta, endpoint=False)
        phi_grid = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
        
        theta_flat, phi_flat = np.meshgrid(theta_grid, phi_grid, indexing='ij')
        theta_flat = theta_flat.flatten()
        phi_flat = phi_flat.flatten()
        
        # Set coordinates
        positions[0, :] = s  # All on same flux surface
        positions[1, :] = theta_flat  # Poloidal angles
        positions[2, :] = phi_flat    # Toroidal angles
        positions[3, :] = velocity    # Uniform velocity
        
        # Pitch angles
        if pitch_angles is not None:
            if len(pitch_angles) == n_particles:
                positions[4, :] = pitch_angles
            else:
                # Cycle through provided pitch angles
                positions[4, :] = np.tile(pitch_angles, (n_particles // len(pitch_angles) + 1))[:n_particles]
        else:
            # Default to zero pitch (perpendicular to field)
            positions[4, :] = 0.0
        
        return batch
    
    def banana_orbits(
        self,
        n_particles: int,
        s: float,
        banana_width: float = 0.1,
        seed: Optional[int] = None
    ) -> ParticleBatch:
        """
        Initialize particles for banana orbit studies.
        
        Args:
            n_particles: Number of particles
            s: Flux surface coordinate
            banana_width: Characteristic banana orbit width
            seed: Random seed
            
        Returns:
            ParticleBatch: Particle batch optimized for banana orbits
        """
        if seed is not None:
            np.random.seed(seed)
        
        batch = ParticleBatch(n_particles)
        positions = batch.positions
        
        # Flux surface with variation for finite banana width
        s_variation = banana_width * (np.random.random(n_particles) - 0.5)
        positions[0, :] = np.clip(s + s_variation, 0.01, 0.99)
        
        # Focus on specific poloidal locations for banana orbits
        # Banana orbits are most pronounced near trapped-passing boundary
        positions[1, :] = np.random.uniform(0, 2*np.pi, n_particles)
        positions[2, :] = np.random.uniform(0, 2*np.pi, n_particles)
        
        # Velocity distribution
        positions[3, :] = 1.0 + 0.2 * (np.random.random(n_particles) - 0.5)
        
        # Pitch angles near trapped-passing boundary
        # Trapped particles have |pitch| < pitch_critical
        critical_pitch = 0.7  # Approximate value
        positions[4, :] = np.random.uniform(-critical_pitch * 0.9, critical_pitch * 0.9, n_particles)
        
        return batch
    
    def energy_pitch_scan(
        self,
        energies: np.ndarray,
        pitch_angles: np.ndarray,
        s: float
    ) -> ParticleBatch:
        """
        Create particles for energy-pitch angle parameter scan.
        
        Args:
            energies: Array of normalized energies
            pitch_angles: Array of pitch angles
            s: Flux surface coordinate
            
        Returns:
            ParticleBatch: Particle batch for parameter scan
        """
        # Create all combinations
        energy_grid, pitch_grid = np.meshgrid(energies, pitch_angles, indexing='ij')
        energy_flat = energy_grid.flatten()
        pitch_flat = pitch_grid.flatten()
        
        n_particles = len(energy_flat)
        batch = ParticleBatch(n_particles)
        positions = batch.positions
        
        # Random spatial distribution on flux surface
        np.random.seed(42)  # Reproducible spatial distribution
        positions[0, :] = s
        positions[1, :] = np.random.uniform(0, 2*np.pi, n_particles)
        positions[2, :] = np.random.uniform(0, 2*np.pi, n_particles)
        
        # Set energies and pitch angles
        positions[3, :] = np.sqrt(energy_flat)  # v = sqrt(2*E/m)
        positions[4, :] = pitch_flat
        
        return batch
    
    def load_from_equilibrium(
        self,
        n_particles: int,
        s: float,
        sampling_method: str = 'uniform'
    ) -> ParticleBatch:
        """
        Load particles using VMEC equilibrium data.
        
        Args:
            n_particles: Number of particles
            s: Flux surface coordinate
            sampling_method: Method for spatial sampling
            
        Returns:
            ParticleBatch: Equilibrium-aware particle batch
        """
        # TODO: Integrate with existing VMEC reader and samplers.f90
        # For now, use uniform sampling with equilibrium-aware constraints
        
        if sampling_method == 'uniform':
            return self.uniform_poloidal(n_particles, s)
        elif sampling_method == 'grid':
            # Choose grid dimensions based on particle count
            n_theta = int(np.sqrt(n_particles))
            n_phi = n_particles // n_theta
            return self.flux_surface_grid(n_theta, n_phi, s)
        else:
            raise ValueError(f"Unknown sampling method: {sampling_method}")
    
    def get_surface_info(self, s: float) -> Dict[str, Any]:
        """
        Get information about flux surface geometry.
        
        Args:
            s: Flux surface coordinate
            
        Returns:
            Dict: Surface geometry information
        """
        # TODO: Extract from VMEC equilibrium
        # For now, return placeholder information
        
        return {
            's': s,
            'vmec_file': str(self.vmec_file),
            'estimated_circumference': 2 * np.pi * np.sqrt(s),  # Approximate
            'aspect_ratio': 3.0,  # Placeholder
            'magnetic_shear': 0.5,  # Placeholder
        }