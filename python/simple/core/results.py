"""
BatchResults class for performance-optimized access to simulation result arrays.

This module provides zero-copy access to existing Fortran result arrays
with vectorized analysis capabilities for large-scale HPC workflows.
"""

import numpy as np
from typing import Optional, Dict, Any, Tuple, NamedTuple
from dataclasses import dataclass
from ..backends.fortran import FortranResultWrapper


@dataclass
class ConfinementStats:
    """Statistics about particle confinement in the simulation"""
    n_total: int                    # Total number of particles
    n_confined: int                 # Number of confined particles
    n_lost: int                     # Number of lost particles
    confined_fraction: float        # Fraction of particles that remain confined
    mean_loss_time: float          # Mean loss time for lost particles
    loss_time_distribution: Tuple[np.ndarray, np.ndarray]  # Histogram of loss times
    
    def __post_init__(self):
        """Validate statistics consistency"""
        if self.n_confined + self.n_lost != self.n_total:
            raise ValueError("Confined + lost particles must equal total")
        
        if not 0 <= self.confined_fraction <= 1:
            raise ValueError("Confined fraction must be between 0 and 1")


class BatchResults:
    """
    Results container wrapping existing Fortran output arrays.
    
    Provides zero-copy access to simulation results with vectorized
    analysis methods optimized for large particle datasets.
    
    Available Arrays:
        - loss_times: Time at which each particle was lost (np.inf if confined)
        - final_positions: Final coordinates of each particle (5, n_particles)
        - trap_parameter: Trapping parameter for each particle
        - perpendicular_invariant: Perpendicular adiabatic invariant
    
    Performance Characteristics:
        - Zero-copy views of existing Fortran arrays
        - Vectorized NumPy operations for analysis
        - Memory-efficient streaming for large results
        - HDF5 export capabilities
    """
    
    def __init__(self, backend_results: FortranResultWrapper):
        """
        Initialize with backend result wrapper.
        
        Args:
            backend_results: Wrapper around Fortran result arrays
        """
        self._backend = backend_results
        self.n_particles = backend_results.n_particles
        self._validate_arrays()
    
    def _validate_arrays(self):
        """Validate result array consistency"""
        try:
            # Check that all arrays have consistent particle count
            arrays_to_check = [
                ('loss_times', self.loss_times),
                ('trap_parameter', self.trap_parameter),
                ('perpendicular_invariant', self.perpendicular_invariant)
            ]
            
            for name, array in arrays_to_check:
                if array.shape[0] != self.n_particles:
                    print(f"Warning: {name} array size {array.shape[0]} != n_particles {self.n_particles}")
            
            # Check final_positions shape
            final_pos = self.final_positions
            expected_shape = (5, self.n_particles)
            if final_pos.shape != expected_shape:
                print(f"Warning: final_positions shape {final_pos.shape} != expected {expected_shape}")
                
        except Exception as e:
            print(f"Array validation warning: {e}")
    
    @property
    def loss_times(self) -> np.ndarray:
        """
        Access to existing times_lost array.
        
        Returns:
            np.ndarray: Shape (n_particles,) with loss times
                        Values are np.inf for confined particles
        """
        return self._backend.get_times_lost_view()
    
    @property
    def final_positions(self) -> np.ndarray:
        """
        Access to existing zend array.
        
        Returns:
            np.ndarray: Shape (5, n_particles) with final coordinates
                - Row 0: final s coordinates
                - Row 1: final theta coordinates
                - Row 2: final phi coordinates
                - Row 3: final v_par coordinates
                - Row 4: final mu coordinates
        """
        return self._backend.get_zend_view()
    
    @property
    def trap_parameter(self) -> np.ndarray:
        """
        Access to existing trap_par array.
        
        Returns:
            np.ndarray: Shape (n_particles,) with trapping parameters
        """
        return self._backend.get_trap_par_view()
    
    @property
    def perpendicular_invariant(self) -> np.ndarray:
        """
        Access to existing perp_inv array.
        
        Returns:
            np.ndarray: Shape (n_particles,) with perpendicular invariants
        """
        return self._backend.get_perp_inv_view()
    
    @property
    def confined_mask(self) -> np.ndarray:
        """
        Boolean mask for confined particles.
        
        Returns:
            np.ndarray: Shape (n_particles,) boolean array
                        True for confined particles (loss_time == inf)
        """
        return self.loss_times == np.inf
    
    @property
    def lost_mask(self) -> np.ndarray:
        """
        Boolean mask for lost particles.
        
        Returns:
            np.ndarray: Shape (n_particles,) boolean array
                        True for lost particles (loss_time < inf)
        """
        return self.loss_times < np.inf
    
    def confinement_statistics(self, time_bins: int = 50) -> ConfinementStats:
        """
        Vectorized analysis of confinement using existing arrays.
        
        Args:
            time_bins: Number of histogram bins for loss time distribution
            
        Returns:
            ConfinementStats: Comprehensive confinement analysis
        """
        lost_mask = self.lost_mask
        confined_mask = self.confined_mask
        
        n_lost = lost_mask.sum()
        n_confined = confined_mask.sum()
        
        # Calculate loss time statistics
        if n_lost > 0:
            lost_times = self.loss_times[lost_mask]
            mean_loss_time = lost_times.mean()
            
            # Create histogram of loss times
            hist_counts, hist_bins = np.histogram(lost_times, bins=time_bins)
            loss_time_distribution = (hist_counts, hist_bins)
        else:
            mean_loss_time = np.inf
            loss_time_distribution = (np.zeros(time_bins), np.linspace(0, 1, time_bins + 1))
        
        return ConfinementStats(
            n_total=self.n_particles,
            n_confined=n_confined,
            n_lost=n_lost,
            confined_fraction=n_confined / self.n_particles,
            mean_loss_time=mean_loss_time,
            loss_time_distribution=loss_time_distribution
        )
    
    def get_confined_particles(self) -> np.ndarray:
        """
        Get final positions for confined particles only.
        
        Returns:
            np.ndarray: Final positions with shape (5, n_confined)
        """
        confined_mask = self.confined_mask
        n_confined = confined_mask.sum()
        
        if n_confined == 0:
            return np.zeros((5, 0))
        
        # Extract confined particle data
        confined_final = self.final_positions[:, confined_mask]
        
        return confined_final
    
    def get_lost_particles(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Get data for lost particles only.
        
        Returns:
            Tuple: (final_positions, loss_times, trap_parameters)
                - final_positions: shape (5, n_lost)
                - loss_times: shape (n_lost,)
                - trap_parameters: shape (n_lost,)
        """
        lost_mask = self.lost_mask
        n_lost = lost_mask.sum()
        
        if n_lost == 0:
            return np.zeros((5, 0)), np.zeros(0), np.zeros(0)
        
        lost_final = self.final_positions[:, lost_mask]
        lost_times = self.loss_times[lost_mask]
        lost_trap_par = self.trap_parameter[lost_mask]
        
        return lost_final, lost_times, lost_trap_par
    
    def analyze_by_surface(self, s_bins: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Analyze confinement statistics by flux surface.
        
        Args:
            s_bins: Array of flux surface bin edges
            
        Returns:
            Dict: Analysis results by surface with keys:
                - 's_centers': Bin centers
                - 'confined_fraction': Fraction confined in each bin
                - 'mean_loss_time': Mean loss time in each bin
                - 'particle_count': Number of particles in each bin
        """
        final_s = self.final_positions[0, :]  # s coordinates
        confined_mask = self.confined_mask
        
        # Digitize particles into s bins
        bin_indices = np.digitize(final_s, s_bins)
        
        n_bins = len(s_bins) - 1
        s_centers = (s_bins[:-1] + s_bins[1:]) / 2
        
        confined_fraction = np.zeros(n_bins)
        mean_loss_time = np.zeros(n_bins)
        particle_count = np.zeros(n_bins, dtype=int)
        
        for i in range(1, n_bins + 1):  # bin_indices start at 1
            mask = bin_indices == i
            count = mask.sum()
            
            if count > 0:
                particle_count[i-1] = count
                confined_fraction[i-1] = confined_mask[mask].sum() / count
                
                lost_in_bin = self.lost_mask[mask]
                if lost_in_bin.any():
                    mean_loss_time[i-1] = self.loss_times[mask][lost_in_bin].mean()
                else:
                    mean_loss_time[i-1] = np.inf
        
        return {
            's_centers': s_centers,
            'confined_fraction': confined_fraction,
            'mean_loss_time': mean_loss_time,
            'particle_count': particle_count
        }
    
    def save_hdf5(self, filename: str, append: bool = False):
        """
        Efficient HDF5 export of all result arrays.
        
        Args:
            filename: Output HDF5 file path
            append: If True, append to existing file
        """
        try:
            import h5py
        except ImportError:
            raise RuntimeError("h5py required for HDF5 export")
        
        mode = 'a' if append else 'w'
        
        with h5py.File(filename, mode) as f:
            # Create group for this batch
            if append:
                group_name = f'batch_{len(f.keys())}'
            else:
                group_name = 'results'
            
            grp = f.create_group(group_name)
            
            # Save main result arrays with compression
            grp.create_dataset('loss_times', data=self.loss_times, 
                             compression='gzip', shuffle=True)
            grp.create_dataset('final_positions', data=self.final_positions,
                             compression='gzip', shuffle=True)
            grp.create_dataset('trap_parameter', data=self.trap_parameter,
                             compression='gzip', shuffle=True)
            grp.create_dataset('perpendicular_invariant', data=self.perpendicular_invariant,
                             compression='gzip', shuffle=True)
            
            # Save metadata
            grp.attrs['n_particles'] = self.n_particles
            grp.attrs['n_confined'] = self.confined_mask.sum()
            grp.attrs['n_lost'] = self.lost_mask.sum()
            grp.attrs['confined_fraction'] = self.confined_mask.sum() / self.n_particles
            
            print(f"Saved {self.n_particles} particle results to {filename}:{group_name}")
    
    def save_numpy(self, filename: str):
        """
        Save results in NumPy .npz format.
        
        Args:
            filename: Output .npz file path
        """
        np.savez_compressed(
            filename,
            loss_times=self.loss_times,
            final_positions=self.final_positions,
            trap_parameter=self.trap_parameter,
            perpendicular_invariant=self.perpendicular_invariant,
            n_particles=self.n_particles
        )
        print(f"Saved results to {filename}")
    
    @classmethod
    def load_numpy(cls, filename: str) -> 'BatchResults':
        """
        Load results from NumPy .npz format.
        
        Args:
            filename: Input .npz file path
            
        Returns:
            BatchResults: Loaded results
        """
        data = np.load(filename)
        
        # Create mock backend for loaded data
        from ..backends.fortran import FortranResultWrapper
        
        class MockBackend(FortranResultWrapper):
            def __init__(self, data_dict):
                self.n_particles = int(data_dict['n_particles'])
                self._data = data_dict
            
            def get_times_lost_view(self):
                return self._data['loss_times']
            
            def get_zend_view(self):
                return self._data['final_positions']
            
            def get_trap_par_view(self):
                return self._data['trap_parameter']
            
            def get_perp_inv_view(self):
                return self._data['perpendicular_invariant']
        
        mock_backend = MockBackend(data)
        return cls(mock_backend)
    
    def summary(self) -> str:
        """
        Generate a summary string of the results.
        
        Returns:
            str: Human-readable summary
        """
        stats = self.confinement_statistics()
        
        summary = f"BatchResults Summary:\n"
        summary += f"  Total particles: {stats.n_total}\n"
        summary += f"  Confined: {stats.n_confined} ({stats.confined_fraction:.1%})\n"
        summary += f"  Lost: {stats.n_lost} ({1-stats.confined_fraction:.1%})\n"
        
        if stats.n_lost > 0:
            summary += f"  Mean loss time: {stats.mean_loss_time:.3f}\n"
        
        return summary
    
    def __repr__(self) -> str:
        """String representation"""
        return f"BatchResults(n_particles={self.n_particles})"
    
    def __str__(self) -> str:
        """Detailed string representation"""
        return self.summary()


# Analysis utilities
def compare_results(results1: BatchResults, results2: BatchResults, 
                   rtol: float = 1e-10, atol: float = 1e-12) -> Dict[str, bool]:
    """
    Compare two BatchResults for golden record validation.
    
    Args:
        results1: First result set
        results2: Second result set  
        rtol: Relative tolerance for floating point comparison
        atol: Absolute tolerance for floating point comparison
        
    Returns:
        Dict: Comparison results with boolean values for each array
    """
    if results1.n_particles != results2.n_particles:
        raise ValueError("Results must have same number of particles")
    
    comparison = {}
    
    # Compare loss times
    comparison['loss_times'] = np.allclose(
        results1.loss_times, results2.loss_times, rtol=rtol, atol=atol
    )
    
    # Compare final positions
    comparison['final_positions'] = np.allclose(
        results1.final_positions, results2.final_positions, rtol=rtol, atol=atol
    )
    
    # Compare trap parameters
    comparison['trap_parameter'] = np.allclose(
        results1.trap_parameter, results2.trap_parameter, rtol=rtol, atol=atol
    )
    
    # Compare perpendicular invariants
    comparison['perpendicular_invariant'] = np.allclose(
        results1.perpendicular_invariant, results2.perpendicular_invariant, rtol=rtol, atol=atol
    )
    
    # Overall match
    comparison['all_match'] = all(comparison.values())
    
    return comparison


def combine_batch_results(results_list: list) -> BatchResults:
    """
    Combine multiple BatchResults into a single result set.
    
    Args:
        results_list: List of BatchResults to combine
        
    Returns:
        BatchResults: Combined results
    """
    if not results_list:
        raise ValueError("Cannot combine empty results list")
    
    total_particles = sum(r.n_particles for r in results_list)
    
    # Combine all arrays
    combined_loss_times = np.concatenate([r.loss_times for r in results_list])
    combined_final_pos = np.concatenate([r.final_positions for r in results_list], axis=1)
    combined_trap_par = np.concatenate([r.trap_parameter for r in results_list])
    combined_perp_inv = np.concatenate([r.perpendicular_invariant for r in results_list])
    
    # Create mock backend for combined data
    class CombinedBackend:
        def __init__(self, n_particles, loss_times, final_pos, trap_par, perp_inv):
            self.n_particles = n_particles
            self._loss_times = loss_times
            self._final_pos = final_pos
            self._trap_par = trap_par
            self._perp_inv = perp_inv
        
        def get_times_lost_view(self):
            return self._loss_times
        
        def get_zend_view(self):
            return self._final_pos
        
        def get_trap_par_view(self):
            return self._trap_par
        
        def get_perp_inv_view(self):
            return self._perp_inv
    
    combined_backend = CombinedBackend(
        total_particles, combined_loss_times, combined_final_pos, 
        combined_trap_par, combined_perp_inv
    )
    
    return BatchResults(combined_backend)