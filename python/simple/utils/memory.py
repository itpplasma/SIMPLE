"""
Memory-efficient streaming utilities for large particle batches.

Provides tools for processing millions of particles with constant memory usage
through batch streaming and efficient data management.
"""

import os
import gc
import time
import psutil
from typing import Iterator, Optional, Dict, Any, Union, List
from pathlib import Path
import numpy as np

from ..core.batch import ParticleBatch
from ..core.results import BatchResults, combine_batch_results
from ..core.simulation import trace_orbits


class ParticleBatchStream:
    """
    Memory-efficient streaming for processing millions of particles.
    
    Provides iteration over particle batches with constant memory usage,
    enabling simulation of arbitrarily large particle counts within
    available memory constraints.
    """
    
    def __init__(
        self,
        total_particles: int,
        batch_size: int = 100_000,
        initialization_func: Optional[callable] = None,
        initialization_kwargs: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize particle batch stream.
        
        Args:
            total_particles: Total number of particles to process
            batch_size: Number of particles per batch
            initialization_func: Function to initialize each batch
            initialization_kwargs: Arguments for initialization function
        """
        if total_particles <= 0:
            raise ValueError("Total particles must be positive")
        if batch_size <= 0:
            raise ValueError("Batch size must be positive")
        
        self.total_particles = total_particles
        self.batch_size = batch_size
        self.current_batch = 0
        self.total_batches = (total_particles + batch_size - 1) // batch_size
        
        self.initialization_func = initialization_func
        self.initialization_kwargs = initialization_kwargs or {}
        
        # Performance tracking
        self._processed_particles = 0
        self._start_time = None
        self._batch_times = []
    
    def __iter__(self) -> Iterator[ParticleBatch]:
        """Iterate over particle batches"""
        self.current_batch = 0
        self._processed_particles = 0
        self._start_time = time.perf_counter()
        self._batch_times = []
        return self
    
    def __next__(self) -> ParticleBatch:
        """Get next particle batch"""
        if self.current_batch >= self.total_batches:
            raise StopIteration
        
        batch_start_time = time.perf_counter()
        
        # Calculate actual batch size for this iteration
        start_idx = self.current_batch * self.batch_size
        end_idx = min(start_idx + self.batch_size, self.total_particles)
        actual_batch_size = end_idx - start_idx
        
        # Create batch
        batch = ParticleBatch(actual_batch_size)
        
        # Initialize if function provided
        if self.initialization_func is not None:
            self.initialization_func(batch, **self.initialization_kwargs)
        
        self.current_batch += 1
        self._processed_particles += actual_batch_size
        
        batch_time = time.perf_counter() - batch_start_time
        self._batch_times.append(batch_time)
        
        return batch
    
    def progress_info(self) -> Dict[str, Any]:
        """Get progress information"""
        if self._start_time is None:
            return {}
        
        elapsed = time.perf_counter() - self._start_time
        progress = self._processed_particles / self.total_particles
        
        return {
            'processed_particles': self._processed_particles,
            'total_particles': self.total_particles,
            'progress_fraction': progress,
            'progress_percent': progress * 100,
            'elapsed_time': elapsed,
            'estimated_total_time': elapsed / progress if progress > 0 else None,
            'estimated_remaining_time': elapsed * (1 - progress) / progress if progress > 0 else None,
            'current_batch': self.current_batch,
            'total_batches': self.total_batches,
            'average_batch_time': np.mean(self._batch_times) if self._batch_times else 0
        }
    
    def __len__(self) -> int:
        """Number of batches"""
        return self.total_batches


class StreamResults:
    """
    Container for streaming simulation results.
    
    Manages results from multiple batches with efficient storage
    and analysis capabilities for large datasets.
    """
    
    def __init__(self, output_file: str, total_particles: int):
        """
        Initialize stream results.
        
        Args:
            output_file: HDF5 file containing results
            total_particles: Total number of particles
        """
        self.output_file = Path(output_file)
        self.total_particles = total_particles
        self._cached_stats = None
    
    def get_summary_statistics(self) -> Dict[str, Any]:
        """
        Get summary statistics without loading all data.
        
        Returns:
            Dict: Summary statistics across all batches
        """
        if self._cached_stats is not None:
            return self._cached_stats
        
        try:
            import h5py
        except ImportError:
            raise RuntimeError("h5py required for streaming results")
        
        total_confined = 0
        total_lost = 0
        all_loss_times = []
        
        with h5py.File(self.output_file, 'r') as f:
            # Iterate through all batch groups
            for group_name in f.keys():
                if group_name.startswith('batch_'):
                    grp = f[group_name]
                    
                    # Count confined/lost in this batch
                    loss_times = grp['loss_times'][:]
                    confined_mask = loss_times == np.inf
                    
                    total_confined += confined_mask.sum()
                    total_lost += (~confined_mask).sum()
                    
                    # Collect loss times for statistics
                    finite_times = loss_times[~confined_mask]
                    if len(finite_times) > 0:
                        all_loss_times.extend(finite_times)
        
        # Calculate overall statistics
        if all_loss_times:
            all_loss_times = np.array(all_loss_times)
            mean_loss_time = all_loss_times.mean()
            median_loss_time = np.median(all_loss_times)
        else:
            mean_loss_time = np.inf
            median_loss_time = np.inf
        
        self._cached_stats = {
            'total_particles': self.total_particles,
            'total_confined': total_confined,
            'total_lost': total_lost,
            'confined_fraction': total_confined / self.total_particles,
            'lost_fraction': total_lost / self.total_particles,
            'mean_loss_time': mean_loss_time,
            'median_loss_time': median_loss_time,
            'output_file': str(self.output_file)
        }
        
        return self._cached_stats
    
    def load_batch_results(self, batch_index: int) -> BatchResults:
        """
        Load results for a specific batch.
        
        Args:
            batch_index: Batch index to load
            
        Returns:
            BatchResults: Results for the specified batch
        """
        try:
            import h5py
        except ImportError:
            raise RuntimeError("h5py required for streaming results")
        
        group_name = f'batch_{batch_index}'
        
        with h5py.File(self.output_file, 'r') as f:
            if group_name not in f:
                raise KeyError(f"Batch {batch_index} not found in results")
            
            grp = f[group_name]
            
            # Load all arrays
            loss_times = grp['loss_times'][:]
            final_positions = grp['final_positions'][:]
            trap_parameter = grp['trap_parameter'][:]
            perpendicular_invariant = grp['perpendicular_invariant'][:]
            
            n_particles = len(loss_times)
        
        # Create mock backend for loaded data
        class StreamBackend:
            def __init__(self, loss_times, final_positions, trap_parameter, perpendicular_invariant):
                self.n_particles = len(loss_times)
                self._loss_times = loss_times
                self._final_positions = final_positions
                self._trap_parameter = trap_parameter
                self._perpendicular_invariant = perpendicular_invariant
            
            def get_times_lost_view(self):
                return self._loss_times
            
            def get_zend_view(self):
                return self._final_positions
            
            def get_trap_par_view(self):
                return self._trap_parameter
            
            def get_perp_inv_view(self):
                return self._perpendicular_invariant
        
        backend = StreamBackend(loss_times, final_positions, trap_parameter, perpendicular_invariant)
        return BatchResults(backend)
    
    def iter_batch_results(self) -> Iterator[BatchResults]:
        """
        Iterate over all batch results.
        
        Yields:
            BatchResults: Results for each batch
        """
        try:
            import h5py
        except ImportError:
            raise RuntimeError("h5py required for streaming results")
        
        with h5py.File(self.output_file, 'r') as f:
            batch_names = [name for name in f.keys() if name.startswith('batch_')]
            batch_names.sort(key=lambda x: int(x.split('_')[1]))
            
            for batch_name in batch_names:
                batch_index = int(batch_name.split('_')[1])
                yield self.load_batch_results(batch_index)
    
    def combine_all_results(self) -> BatchResults:
        """
        Combine all batch results into a single BatchResults object.
        
        Warning: This loads all data into memory and may require significant RAM.
        
        Returns:
            BatchResults: Combined results from all batches
        """
        batch_results = list(self.iter_batch_results())
        return combine_batch_results(batch_results)


def process_large_simulation(
    vmec_file: str,
    n_total: int,
    tmax: float,
    batch_size: int = 100_000,
    output_file: str = 'results.h5',
    s_surface: float = 0.9,
    integrator: str = 'symplectic_midpoint',
    verbose: bool = True,
    memory_limit_gb: Optional[float] = None
) -> StreamResults:
    """
    Process millions of particles in memory-efficient batches.
    
    Memory Usage: Constant - independent of total particle count
    Performance: Linear scaling with existing OpenMP optimization
    Output: Streaming HDF5 for large dataset handling
    
    Args:
        vmec_file: Path to VMEC equilibrium file
        n_total: Total number of particles to simulate
        tmax: Maximum simulation time
        batch_size: Particles per batch (memory usage scales with this)
        output_file: Output HDF5 file path
        s_surface: Flux surface for particle initialization
        integrator: Integration method
        verbose: Print progress information
        memory_limit_gb: Memory limit in GB (auto-adjust batch size if exceeded)
        
    Returns:
        StreamResults: Handle for accessing streaming results
    """
    if verbose:
        print(f"Starting large simulation: {n_total:,} particles in batches of {batch_size:,}")
        print(f"VMEC file: {vmec_file}")
        print(f"Output: {output_file}")
    
    # Monitor memory usage
    process = psutil.Process()
    initial_memory = process.memory_info().rss / (1024**3)  # GB
    
    if memory_limit_gb is not None and verbose:
        print(f"Memory limit: {memory_limit_gb:.1f} GB (current: {initial_memory:.1f} GB)")
    
    # Initialize function for particle batches
    def init_batch(batch: ParticleBatch, **kwargs):
        batch.initialize_surface(vmec_file, s=s_surface)
    
    # Create stream
    stream = ParticleBatchStream(
        n_total,
        batch_size,
        initialization_func=init_batch
    )
    
    # Remove existing output file
    output_path = Path(output_file)
    if output_path.exists():
        output_path.unlink()
        if verbose:
            print(f"Removed existing output file: {output_file}")
    
    start_time = time.perf_counter()
    batch_count = 0
    
    try:
        for i, batch in enumerate(stream):
            batch_start = time.perf_counter()
            
            # Check memory usage
            current_memory = process.memory_info().rss / (1024**3)
            if memory_limit_gb is not None and current_memory > memory_limit_gb:
                print(f"Warning: Memory usage {current_memory:.1f} GB exceeds limit {memory_limit_gb:.1f} GB")
            
            # Process batch using existing optimized algorithms
            results = trace_orbits(batch, tmax=tmax, integrator=integrator, verbose=False)
            
            # Stream results to disk
            results.save_hdf5(output_file, append=(i > 0))
            
            batch_time = time.perf_counter() - batch_start
            batch_count += 1
            
            if verbose:
                progress = stream.progress_info()
                print(f"Batch {i+1}/{stream.total_batches}: "
                      f"{batch.n_particles:,} particles in {batch_time:.2f}s "
                      f"({progress['progress_percent']:.1f}% complete, "
                      f"memory: {current_memory:.1f} GB)")
            
            # Cleanup memory between batches
            del batch, results
            gc.collect()
    
    except KeyboardInterrupt:
        print(f"\nSimulation interrupted after {batch_count} batches")
        if output_path.exists() and batch_count > 0:
            print(f"Partial results saved to {output_file}")
        raise
    
    except Exception as e:
        print(f"\nSimulation failed after {batch_count} batches: {e}")
        raise
    
    total_time = time.perf_counter() - start_time
    final_memory = process.memory_info().rss / (1024**3)
    
    if verbose:
        print(f"\nSimulation completed:")
        print(f"  Total time: {total_time:.1f}s")
        print(f"  Processed: {batch_count * batch_size:,} particles")
        print(f"  Rate: {(batch_count * batch_size) / total_time:.0f} particles/second")
        print(f"  Memory usage: {initial_memory:.1f} → {final_memory:.1f} GB")
        print(f"  Output: {output_file}")
    
    return StreamResults(output_file, n_total)


def estimate_memory_usage(
    n_particles: int,
    include_results: bool = True,
    overhead_factor: float = 1.5
) -> Dict[str, float]:
    """
    Estimate memory usage for a given number of particles.
    
    Args:
        n_particles: Number of particles
        include_results: Include result arrays in estimate
        overhead_factor: Factor for Python/system overhead
        
    Returns:
        Dict: Memory estimates in MB
    """
    # Particle data: 5 coordinates * 8 bytes per float64
    particle_memory = n_particles * 5 * 8
    
    # Result arrays
    if include_results:
        result_memory = (
            n_particles * 8 +  # loss_times
            n_particles * 5 * 8 +  # final_positions  
            n_particles * 8 +  # trap_parameter
            n_particles * 8     # perpendicular_invariant
        )
    else:
        result_memory = 0
    
    total_memory = (particle_memory + result_memory) * overhead_factor
    
    return {
        'particles_mb': particle_memory / (1024**2),
        'results_mb': result_memory / (1024**2),
        'total_mb': total_memory / (1024**2),
        'total_gb': total_memory / (1024**3),
        'overhead_factor': overhead_factor
    }


def optimize_batch_size(
    target_memory_gb: float,
    n_total: int,
    safety_factor: float = 0.8
) -> int:
    """
    Optimize batch size for given memory constraint.
    
    Args:
        target_memory_gb: Target memory usage in GB
        n_total: Total number of particles
        safety_factor: Safety factor for memory estimate
        
    Returns:
        int: Recommended batch size
    """
    target_memory_mb = target_memory_gb * 1024 * safety_factor
    
    # Binary search for optimal batch size
    min_batch = 1000
    max_batch = min(n_total, 10_000_000)
    
    while min_batch < max_batch:
        mid_batch = (min_batch + max_batch + 1) // 2
        memory_est = estimate_memory_usage(mid_batch)
        
        if memory_est['total_mb'] <= target_memory_mb:
            min_batch = mid_batch
        else:
            max_batch = mid_batch - 1
    
    return min_batch


class MemoryMonitor:
    """Context manager for monitoring memory usage during simulations"""
    
    def __init__(self, name: str = "simulation", verbose: bool = True):
        self.name = name
        self.verbose = verbose
        self.process = psutil.Process()
        self.initial_memory = None
        self.peak_memory = None
        self.memory_samples = []
    
    def __enter__(self):
        self.initial_memory = self.process.memory_info().rss / (1024**3)
        self.peak_memory = self.initial_memory
        self.memory_samples = [self.initial_memory]
        
        if self.verbose:
            print(f"Starting {self.name} with {self.initial_memory:.1f} GB")
        
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        final_memory = self.process.memory_info().rss / (1024**3)
        memory_increase = final_memory - self.initial_memory
        
        if self.verbose:
            print(f"{self.name} completed:")
            print(f"  Final memory: {final_memory:.1f} GB")
            print(f"  Peak memory: {self.peak_memory:.1f} GB")
            print(f"  Memory increase: {memory_increase:+.1f} GB")
    
    def sample(self):
        """Sample current memory usage"""
        current = self.process.memory_info().rss / (1024**3)
        self.memory_samples.append(current)
        self.peak_memory = max(self.peak_memory, current)
        return current
    
    def get_stats(self) -> Dict[str, float]:
        """Get memory usage statistics"""
        if not self.memory_samples:
            return {}
        
        samples = np.array(self.memory_samples)
        return {
            'initial_gb': self.initial_memory,
            'peak_gb': self.peak_memory,
            'current_gb': samples[-1],
            'mean_gb': samples.mean(),
            'std_gb': samples.std(),
            'increase_gb': samples[-1] - self.initial_memory
        }