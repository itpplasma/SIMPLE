"""
Memory utilities for large-scale processing using existing Fortran batch capabilities.

Pure interface layer - all batch processing delegated to existing Fortran implementations.
No independent streaming logic - leverages existing batch processing patterns.
"""

import os
import gc
import psutil
from typing import Optional, Dict, Any
import numpy as np

from ..core.batch import ParticleBatch
from ..core.results import BatchResults
from ..core.simulation import trace_orbits


def estimate_memory_usage(n_particles: int) -> Dict[str, float]:
    """
    Estimate memory usage for particle processing.
    
    Args:
        n_particles: Number of particles
        
    Returns:
        Dict: Memory estimates in MB
    """
    # Based on existing Fortran array sizes
    # zstart(5, n_particles) + results arrays
    
    particle_data_mb = (5 * n_particles * 8) / (1024 * 1024)  # 8 bytes per double
    result_arrays_mb = (3 * n_particles * 8) / (1024 * 1024)  # times_lost, trap_par, perp_inv
    overhead_mb = particle_data_mb * 0.2  # Estimated overhead
    
    total_mb = particle_data_mb + result_arrays_mb + overhead_mb
    
    return {
        'particle_data_mb': particle_data_mb,
        'result_arrays_mb': result_arrays_mb,
        'overhead_mb': overhead_mb,
        'total_mb': total_mb
    }


def optimize_batch_size(target_memory_mb: float = 1000.0) -> int:
    """
    Calculate optimal batch size for available memory.
    
    Args:
        target_memory_mb: Target memory usage in MB
        
    Returns:
        int: Recommended batch size
    """
    # Binary search for optimal batch size
    low, high = 1000, 1_000_000
    
    while low < high:
        mid = (low + high + 1) // 2
        mem_usage = estimate_memory_usage(mid)
        
        if mem_usage['total_mb'] <= target_memory_mb:
            low = mid
        else:
            high = mid - 1
    
    return low


def get_available_memory_mb() -> float:
    """
    Get available system memory in MB.
    
    Returns:
        float: Available memory in MB
    """
    try:
        return psutil.virtual_memory().available / (1024 * 1024)
    except:
        # Fallback estimate if psutil unavailable
        return 4000.0  # Conservative 4GB estimate


def process_large_simulation(
    total_particles: int,
    vmec_file: str,
    tmax: float,
    max_memory_mb: float = 2000.0,
    **simulation_kwargs
) -> BatchResults:
    """
    Process large simulation using existing Fortran batch processing.
    
    This function coordinates multiple calls to existing trace_orbits()
    without implementing independent streaming logic.
    
    Args:
        total_particles: Total number of particles
        vmec_file: VMEC equilibrium file
        tmax: Simulation time
        max_memory_mb: Maximum memory to use
        **simulation_kwargs: Additional arguments for trace_orbits()
        
    Returns:
        BatchResults: Combined results from all batches
    """
    # Calculate optimal batch size using existing memory estimation
    batch_size = optimize_batch_size(max_memory_mb)
    
    if total_particles <= batch_size:
        # Single batch - delegate directly to existing implementation
        from ..samplers import SurfaceSampler
        sampler = SurfaceSampler(vmec_file)
        particle_data = sampler.sample_surface_fieldline(total_particles)
        particles = ParticleBatch.from_fortran_arrays(particle_data)
        
        return trace_orbits(particles, tmax=tmax, **simulation_kwargs)
    
    # Multiple batches - coordinate calls to existing trace_orbits()
    print(f"Processing {total_particles} particles in batches of {batch_size}")
    
    all_results = []
    processed = 0
    
    while processed < total_particles:
        current_batch_size = min(batch_size, total_particles - processed)
        
        # Use existing sampling for each batch
        from ..samplers import SurfaceSampler
        sampler = SurfaceSampler(vmec_file)
        particle_data = sampler.sample_surface_fieldline(current_batch_size)
        particles = ParticleBatch.from_fortran_arrays(particle_data)
        
        # Delegate to existing trace_orbits implementation
        batch_results = trace_orbits(particles, tmax=tmax, **simulation_kwargs)
        all_results.append(batch_results)
        
        processed += current_batch_size
        print(f"Completed batch: {processed}/{total_particles} particles")
        
        # Force garbage collection between batches
        gc.collect()
    
    # Combine results using existing array operations
    from ..core.results import combine_batch_results
    return combine_batch_results(all_results)


class MemoryMonitor:
    """Simple memory monitoring for batch processing"""
    
    def __init__(self):
        self.initial_memory = self._get_memory_usage()
        self.peak_memory = self.initial_memory
    
    def _get_memory_usage(self) -> float:
        """Get current memory usage in MB"""
        try:
            process = psutil.Process(os.getpid())
            return process.memory_info().rss / (1024 * 1024)
        except:
            return 0.0
    
    def update(self):
        """Update peak memory tracking"""
        current = self._get_memory_usage()
        self.peak_memory = max(self.peak_memory, current)
    
    def get_usage(self) -> Dict[str, float]:
        """Get memory usage statistics"""
        current = self._get_memory_usage()
        return {
            'initial_mb': self.initial_memory,
            'current_mb': current,
            'peak_mb': self.peak_memory,
            'delta_mb': current - self.initial_memory
        }


# Convenience function for common large-scale processing
def run_million_particle_simulation(
    vmec_file: str,
    n_million: float = 1.0,
    tmax: float = 1000.0,
    **kwargs
) -> BatchResults:
    """
    Run simulation with millions of particles using existing batch processing.
    
    Args:
        vmec_file: VMEC equilibrium file
        n_million: Number of millions of particles (e.g., 1.5 = 1.5M particles)
        tmax: Simulation time
        **kwargs: Additional simulation parameters
        
    Returns:
        BatchResults: Combined simulation results
    """
    total_particles = int(n_million * 1_000_000)
    
    # Use existing large-scale processing coordination
    return process_large_simulation(
        total_particles=total_particles,
        vmec_file=vmec_file,
        tmax=tmax,
        **kwargs
    )