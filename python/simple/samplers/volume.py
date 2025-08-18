"""
Volume sampling for particle initialization between flux surfaces.

Provides various volume sampling methods for comprehensive plasma
simulations using existing SIMPLE volume sampling algorithms.
"""

import numpy as np
from typing import Optional, Tuple, Dict, Any
from pathlib import Path

from ..core.batch import ParticleBatch


class VolumeSampler:
    """
    Volume sampling for particle initialization between flux surfaces.
    
    Provides methods for distributing particles throughout a volume
    between specified flux surfaces with various distribution patterns.
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
    
    def uniform_volume(
        self,
        n_particles: int,
        s_range: Tuple[float, float],
        velocity_range: Tuple[float, float] = (0.5, 1.5),
        pitch_range: Tuple[float, float] = (-1.0, 1.0),
        seed: Optional[int] = None
    ) -> ParticleBatch:
        """
        Uniform distribution throughout volume between flux surfaces.
        
        Args:
            n_particles: Number of particles to generate
            s_range: (s_min, s_max) flux surface range
            velocity_range: (v_min, v_max) velocity range
            pitch_range: (pitch_min, pitch_max) pitch angle range
            seed: Random seed for reproducibility
            
        Returns:
            ParticleBatch: Volume-initialized particle batch
        """
        s_min, s_max = s_range
        if not (0 < s_min < s_max < 1):
            raise ValueError("Must have 0 < s_min < s_max < 1")
        
        if seed is not None:
            np.random.seed(seed)
        
        batch = ParticleBatch(n_particles)
        positions = batch.positions
        
        # Uniform flux coordinate distribution
        # For true volume uniformity, should weight by volume element
        positions[0, :] = np.random.uniform(s_min, s_max, n_particles)
        
        # Uniform angular distributions
        positions[1, :] = np.random.uniform(0, 2*np.pi, n_particles)  # theta
        positions[2, :] = np.random.uniform(0, 2*np.pi, n_particles)  # phi
        
        # Velocity distribution
        v_min, v_max = velocity_range
        positions[3, :] = np.random.uniform(v_min, v_max, n_particles)
        
        # Pitch angle distribution
        pitch_min, pitch_max = pitch_range
        positions[4, :] = np.random.uniform(pitch_min, pitch_max, n_particles)
        
        return batch
    
    def volume_weighted(
        self,
        n_particles: int,
        s_range: Tuple[float, float],
        seed: Optional[int] = None
    ) -> ParticleBatch:
        """
        Volume-weighted distribution accounting for flux surface expansion.
        
        Args:
            n_particles: Number of particles
            s_range: (s_min, s_max) flux surface range
            seed: Random seed
            
        Returns:
            ParticleBatch: Volume-weighted particle batch
        """
        s_min, s_max = s_range
        if not (0 < s_min < s_max < 1):
            raise ValueError("Must have 0 < s_min < s_max < 1")
        
        if seed is not None:
            np.random.seed(seed)
        
        batch = ParticleBatch(n_particles)
        positions = batch.positions
        
        # Volume-weighted s distribution
        # For toroidal geometry: dV ~ s ds (approximate)
        # Use inverse transform sampling: s = sqrt(s_min^2 + u*(s_max^2 - s_min^2))
        u = np.random.random(n_particles)
        s_squared = s_min**2 + u * (s_max**2 - s_min**2)
        positions[0, :] = np.sqrt(s_squared)
        
        # Uniform angular distributions
        positions[1, :] = np.random.uniform(0, 2*np.pi, n_particles)
        positions[2, :] = np.random.uniform(0, 2*np.pi, n_particles)
        
        # Default velocity and pitch
        positions[3, :] = 1.0 + 0.1 * (np.random.random(n_particles) - 0.5)
        positions[4, :] = np.random.uniform(-1.0, 1.0, n_particles)
        
        return batch
    
    def radial_profile(
        self,
        n_particles: int,
        s_range: Tuple[float, float],
        profile_type: str = 'parabolic',
        profile_params: Optional[Dict[str, float]] = None,
        seed: Optional[int] = None
    ) -> ParticleBatch:
        """
        Distribution following specified radial density profile.
        
        Args:
            n_particles: Number of particles
            s_range: (s_min, s_max) flux surface range
            profile_type: 'parabolic', 'exponential', 'linear', or 'custom'
            profile_params: Parameters for profile shape
            seed: Random seed
            
        Returns:
            ParticleBatch: Profile-distributed particle batch
        """
        s_min, s_max = s_range
        if seed is not None:
            np.random.seed(seed)
        
        if profile_params is None:
            profile_params = {}
        
        # Generate s coordinates according to profile
        if profile_type == 'parabolic':
            # n(s) = n0 * (1 - s^alpha)
            alpha = profile_params.get('alpha', 2.0)
            s_coords = self._sample_parabolic_profile(n_particles, s_min, s_max, alpha)
        elif profile_type == 'exponential':
            # n(s) = n0 * exp(-s/lambda)
            lambda_s = profile_params.get('lambda', 0.3)
            s_coords = self._sample_exponential_profile(n_particles, s_min, s_max, lambda_s)
        elif profile_type == 'linear':
            # n(s) = n0 * (1 - s/s_edge)
            s_coords = self._sample_linear_profile(n_particles, s_min, s_max)
        else:
            raise ValueError(f"Unknown profile type: {profile_type}")
        
        batch = ParticleBatch(n_particles)
        positions = batch.positions
        
        positions[0, :] = s_coords
        positions[1, :] = np.random.uniform(0, 2*np.pi, n_particles)
        positions[2, :] = np.random.uniform(0, 2*np.pi, n_particles)
        positions[3, :] = 1.0 + 0.1 * (np.random.random(n_particles) - 0.5)
        positions[4, :] = np.random.uniform(-1.0, 1.0, n_particles)
        
        return batch
    
    def _sample_parabolic_profile(self, n_particles: int, s_min: float, s_max: float, alpha: float) -> np.ndarray:
        """Sample from parabolic density profile using rejection sampling"""
        # For simplicity, use inverse transform for alpha=2 case
        if abs(alpha - 2.0) < 1e-6:
            u = np.random.random(n_particles)
            # For n(s) = 1 - s^2, CDF requires numerical integration
            # Approximate with simple power law
            s_coords = s_min + (s_max - s_min) * (1 - u**(1/3))
        else:
            # General case - use rejection sampling
            s_coords = []
            while len(s_coords) < n_particles:
                s_test = np.random.uniform(s_min, s_max, n_particles * 2)
                prob_test = (1 - s_test**alpha) / (1 - s_min**alpha)
                accept = np.random.random(len(s_test)) < prob_test
                s_coords.extend(s_test[accept])
            s_coords = np.array(s_coords[:n_particles])
        
        return s_coords
    
    def _sample_exponential_profile(self, n_particles: int, s_min: float, s_max: float, lambda_s: float) -> np.ndarray:
        """Sample from exponential density profile"""
        # Use inverse transform sampling for exponential
        u = np.random.random(n_particles)
        # Truncated exponential between s_min and s_max
        exp_min = np.exp(-s_min / lambda_s)
        exp_max = np.exp(-s_max / lambda_s)
        exp_vals = exp_min + u * (exp_max - exp_min)
        s_coords = -lambda_s * np.log(exp_vals)
        
        return np.clip(s_coords, s_min, s_max)
    
    def _sample_linear_profile(self, n_particles: int, s_min: float, s_max: float) -> np.ndarray:
        """Sample from linear density profile"""
        # For n(s) = 1 - s, use inverse transform
        u = np.random.random(n_particles)
        # Analytical inverse for linear profile
        s_coords = s_min + (s_max - s_min) * np.sqrt(u)
        
        return s_coords
    
    def energy_distribution(
        self,
        n_particles: int,
        s_range: Tuple[float, float],
        energy_profile: str = 'maxwellian',
        temperature: float = 1.0,
        seed: Optional[int] = None
    ) -> ParticleBatch:
        """
        Volume sampling with specified energy distribution.
        
        Args:
            n_particles: Number of particles
            s_range: (s_min, s_max) flux surface range
            energy_profile: 'maxwellian', 'monoenergetic', or 'power_law'
            temperature: Characteristic temperature for distribution
            seed: Random seed
            
        Returns:
            ParticleBatch: Energy-distributed particle batch
        """
        if seed is not None:
            np.random.seed(seed)
        
        batch = ParticleBatch(n_particles)
        positions = batch.positions
        
        # Uniform spatial distribution
        s_min, s_max = s_range
        positions[0, :] = np.random.uniform(s_min, s_max, n_particles)
        positions[1, :] = np.random.uniform(0, 2*np.pi, n_particles)
        positions[2, :] = np.random.uniform(0, 2*np.pi, n_particles)
        
        # Energy distribution
        if energy_profile == 'maxwellian':
            # Maxwell-Boltzmann velocity distribution
            # v = sqrt(2*E/m) where E ~ Gamma(3/2, kT)
            energies = np.random.gamma(1.5, temperature, n_particles)
            velocities = np.sqrt(2 * energies)
        elif energy_profile == 'monoenergetic':
            # All particles have same energy
            velocities = np.full(n_particles, np.sqrt(2 * temperature))
        elif energy_profile == 'power_law':
            # Power law energy distribution E^(-alpha)
            alpha = 3.0  # Typical for fusion plasmas
            u = np.random.random(n_particles)
            energies = temperature * ((1 - u) ** (-1/(alpha - 1)))
            velocities = np.sqrt(2 * energies)
        else:
            raise ValueError(f"Unknown energy profile: {energy_profile}")
        
        positions[3, :] = velocities
        
        # Isotropic pitch angle distribution
        positions[4, :] = np.random.uniform(-1.0, 1.0, n_particles)
        
        return batch
    
    def layered_sampling(
        self,
        n_particles: int,
        s_layers: np.ndarray,
        particles_per_layer: Optional[np.ndarray] = None,
        seed: Optional[int] = None
    ) -> ParticleBatch:
        """
        Sample particles in discrete radial layers.
        
        Args:
            n_particles: Total number of particles
            s_layers: Array of flux surface coordinates for layers
            particles_per_layer: Particles in each layer (uniform if None)
            seed: Random seed
            
        Returns:
            ParticleBatch: Layer-sampled particle batch
        """
        if seed is not None:
            np.random.seed(seed)
        
        n_layers = len(s_layers)
        
        if particles_per_layer is None:
            # Distribute evenly across layers
            base_count = n_particles // n_layers
            particles_per_layer = np.full(n_layers, base_count)
            # Distribute remainder
            remainder = n_particles % n_layers
            particles_per_layer[:remainder] += 1
        else:
            if len(particles_per_layer) != n_layers:
                raise ValueError("particles_per_layer must have same length as s_layers")
            if sum(particles_per_layer) != n_particles:
                raise ValueError("Sum of particles_per_layer must equal n_particles")
        
        all_particles = []
        
        for i, (s_layer, n_layer) in enumerate(zip(s_layers, particles_per_layer)):
            if n_layer == 0:
                continue
            
            # Create particles for this layer
            layer_positions = np.zeros((5, n_layer))
            
            # Small variation around layer s coordinate
            layer_positions[0, :] = s_layer + 0.01 * (np.random.random(n_layer) - 0.5)
            layer_positions[1, :] = np.random.uniform(0, 2*np.pi, n_layer)
            layer_positions[2, :] = np.random.uniform(0, 2*np.pi, n_layer)
            layer_positions[3, :] = 1.0 + 0.1 * (np.random.random(n_layer) - 0.5)
            layer_positions[4, :] = np.random.uniform(-1.0, 1.0, n_layer)
            
            all_particles.append(layer_positions)
        
        # Combine all layers
        combined_positions = np.concatenate(all_particles, axis=1)
        
        batch = ParticleBatch(n_particles)
        batch.positions[:] = combined_positions
        
        return batch
    
    def get_volume_info(self, s_range: Tuple[float, float]) -> Dict[str, Any]:
        """
        Get information about volume between flux surfaces.
        
        Args:
            s_range: (s_min, s_max) flux surface range
            
        Returns:
            Dict: Volume geometry information
        """
        s_min, s_max = s_range
        
        # TODO: Calculate from VMEC equilibrium
        # For now, return approximate values
        
        volume_approx = np.pi**2 * (s_max**2 - s_min**2)  # Approximate toroidal volume
        
        return {
            's_range': s_range,
            'vmec_file': str(self.vmec_file),
            'estimated_volume': volume_approx,
            'aspect_ratio': 3.0,  # Placeholder
            'average_s': (s_min + s_max) / 2,
            'volume_fraction': s_max**2 - s_min**2,  # Fraction of total volume
        }