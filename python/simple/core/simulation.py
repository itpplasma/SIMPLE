"""
Pure interface to existing Fortran simulation functionality.

This module provides zero-copy access to existing simple_main.f90 execution
without implementing any computational logic.
"""

import os
import time
from typing import Dict, Any, Optional, Union
from pathlib import Path
import numpy as np

from .batch import ParticleBatch
from .results import BatchResults
from ..backends.fortran import get_backend


# Integrator mapping to existing Fortran implementation
_INTEGRATOR_MAP = {
    'symplectic': 1,
    'symplectic_euler': 1,
    'symplectic_midpoint': 2,
    'symplectic_gauss': 3,
    'symplectic_lobatto': 4,
    'quasi_symplectic': 5,
    'rk45': 6
}


def trace_orbits(
    particles: ParticleBatch,
    tmax: float,
    integrator: str = 'symplectic_midpoint',
    openmp_threads: Optional[int] = None,
    config: Optional[Dict[str, Any]] = None,
    verbose: bool = False
) -> BatchResults:
    """
    Delegate to existing simple_main.f90 execution via f90wrap.
    
    Pure interface - no independent computation logic.
    All orbital integration delegated to proven Fortran implementation.
    
    Args:
        particles: ParticleBatch with initialized particle data
        tmax: Maximum simulation time
        integrator: Integration method name
        openmp_threads: Number of OpenMP threads (None = auto)
        config: Pre-created configuration (overrides other parameters)
        verbose: Print execution information
        
    Returns:
        BatchResults: Zero-copy wrapper around Fortran result arrays
    """
    start_time = time.perf_counter()
    
    # Configure OpenMP threads if specified
    original_omp_threads = None
    if openmp_threads is not None:
        original_omp_threads = os.environ.get('OMP_NUM_THREADS')
        os.environ['OMP_NUM_THREADS'] = str(openmp_threads)
        if verbose:
            print(f"Set OMP_NUM_THREADS = {openmp_threads}")
    
    try:
        # Create configuration for existing Fortran execution
        if config is None:
            config = {
                'tmax': tmax,
                'integmode': _INTEGRATOR_MAP.get(integrator, 2),
                'ntestpart': particles.n_particles,
            }
        else:
            config = config.copy()
            config['ntestpart'] = particles.n_particles
        
        if verbose:
            print(f"Delegating to Fortran with config: {config}")
            print(f"Particle batch: {particles.n_particles} particles")
        
        # Get backend and delegate to existing simple_main execution
        backend = get_backend()
        
        # Execute using existing optimized Fortran implementation
        # This calls simple_main.f90 run() subroutine via f90wrap
        result_wrapper = backend.run_simulation(particles._arrays, config)
        
        # Create results wrapper around existing Fortran arrays
        results = BatchResults(result_wrapper)
        
        total_time = time.perf_counter() - start_time
        
        if verbose:
            print(f"Fortran execution completed in {total_time:.3f}s")
        
        return results
        
    finally:
        # Restore original OpenMP thread setting
        if original_omp_threads is not None:
            os.environ['OMP_NUM_THREADS'] = original_omp_threads
        elif openmp_threads is not None:
            os.environ.pop('OMP_NUM_THREADS', None)


def create_configuration(
    vmec_file: str,
    tmax: float,
    integrator: str = 'symplectic_midpoint',
    dtau: Optional[float] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Create configuration for existing params.f90 namelist.
    
    Args:
        vmec_file: Path to VMEC equilibrium file
        tmax: Maximum simulation time
        integrator: Integration method name
        dtau: Time step (uses Fortran defaults if None)
        **kwargs: Additional parameters for existing namelist
        
    Returns:
        Dict: Configuration dictionary for Fortran execution
    """
    if integrator not in _INTEGRATOR_MAP:
        raise ValueError(f"Unknown integrator '{integrator}'. "
                        f"Available: {list(_INTEGRATOR_MAP.keys())}")
    
    vmec_path = Path(vmec_file)
    if not vmec_path.exists():
        raise FileNotFoundError(f"VMEC file not found: {vmec_file}")
    
    config = {
        'vmec_file': str(vmec_path.absolute()),
        'tmax': float(tmax),
        'integmode': _INTEGRATOR_MAP[integrator],
        **kwargs
    }
    
    # Add time step if specified (otherwise use Fortran defaults)
    if dtau is not None:
        config['dtau'] = float(dtau)
    
    return config


def quick_simulation(
    n_particles: int,
    vmec_file: str = "wout.nc",
    s_surface: float = 0.9,
    tmax: float = 1000.0,
    integrator: str = 'symplectic_midpoint'
) -> BatchResults:
    """
    Quick simulation using existing sampling and execution.
    
    Args:
        n_particles: Number of particles
        vmec_file: VMEC equilibrium file
        s_surface: Flux surface for particle initialization
        tmax: Maximum simulation time
        integrator: Integration method
        
    Returns:
        BatchResults: Simulation results
    """
    # Use existing sampling to initialize particles
    from ..samplers import SurfaceSampler
    
    sampler = SurfaceSampler(vmec_file)
    particle_data = sampler.sample_surface_fieldline(n_particles)
    
    # Create batch from existing Fortran arrays
    particles = ParticleBatch.from_fortran_arrays(particle_data)
    
    # Delegate to existing execution
    results = trace_orbits(particles, tmax=tmax, integrator=integrator, verbose=True)
    
    return results


def validate_golden_record(
    reference_file: str,
    test_particles: ParticleBatch,
    config: Dict[str, Any],
    rtol: float = 1e-10,
    atol: float = 1e-12
) -> Dict[str, bool]:
    """
    Validate against golden record using existing comparison.
    
    Args:
        reference_file: Path to reference results file
        test_particles: Particles to simulate
        config: Simulation configuration
        rtol: Relative tolerance
        atol: Absolute tolerance
        
    Returns:
        Dict: Validation results
    """
    # Load reference using existing format
    reference = BatchResults.load_numpy(reference_file)
    
    # Run test using existing execution
    integrator_name = next(k for k, v in _INTEGRATOR_MAP.items() if v == config['integmode'])
    test_results = trace_orbits(
        test_particles,
        tmax=config['tmax'],
        integrator=integrator_name,
        config=config
    )
    
    # Compare using existing comparison logic
    comparison = reference.compare_with(test_results, rtol=rtol, atol=atol)
    
    return comparison


def benchmark_performance(
    n_particles: int,
    tmax: float = 100.0,
    integrator: str = 'symplectic_midpoint',
    n_runs: int = 3,
    openmp_threads: Optional[int] = None
) -> Dict[str, float]:
    """
    Benchmark interface overhead vs pure Fortran execution.
    
    Args:
        n_particles: Number of particles for benchmark
        tmax: Simulation time
        integrator: Integration method
        n_runs: Number of benchmark runs
        openmp_threads: OpenMP thread count
        
    Returns:
        Dict: Performance metrics
    """
    print(f"Benchmarking interface overhead with {n_particles} particles...")
    
    # Use existing sampling for initialization
    from ..samplers import SurfaceSampler
    sampler = SurfaceSampler("wout.nc")
    particle_data = sampler.sample_surface_fieldline(n_particles)
    batch = ParticleBatch.from_fortran_arrays(particle_data)
    
    times = []
    
    for run in range(n_runs):
        start_time = time.perf_counter()
        
        try:
            # Delegate to existing Fortran execution
            results = trace_orbits(
                batch, 
                tmax=tmax, 
                integrator=integrator,
                openmp_threads=openmp_threads
            )
            
            end_time = time.perf_counter()
            times.append(end_time - start_time)
            
            print(f"Run {run+1}: {times[-1]:.3f}s")
            
        except Exception as e:
            print(f"Run {run+1} failed: {e}")
            continue
    
    if not times:
        raise RuntimeError("All benchmark runs failed")
    
    metrics = {
        'mean_time': np.mean(times),
        'std_time': np.std(times),
        'min_time': np.min(times),
        'max_time': np.max(times),
        'particles_per_second': n_particles / np.mean(times),
    }
    
    print(f"Interface benchmark results:")
    print(f"  Mean time: {metrics['mean_time']:.3f} ± {metrics['std_time']:.3f}s")
    print(f"  Particles/sec: {metrics['particles_per_second']:.0f}")
    
    return metrics


def parameter_sweep(
    base_config: Dict[str, Any],
    parameter_name: str,
    parameter_values: list,
    particles: ParticleBatch
) -> Dict[Any, BatchResults]:
    """
    Parameter sweep using existing execution for each configuration.
    
    Args:
        base_config: Base simulation configuration
        parameter_name: Parameter to sweep
        parameter_values: List of parameter values to test
        particles: Particle batch to use for all runs
        
    Returns:
        Dict: Results for each parameter value
    """
    results = {}
    
    for value in parameter_values:
        config = base_config.copy()
        config[parameter_name] = value
        
        print(f"Running with {parameter_name} = {value}")
        
        # Use fresh copy for each run
        test_particles = particles.copy()
        
        # Delegate to existing execution
        result = trace_orbits(
            test_particles,
            tmax=config['tmax'],
            config=config
        )
        
        results[value] = result
        
        # Summary from existing result functionality
        stats = result.confinement_statistics()
        print(f"  Confined fraction: {stats.confined_fraction:.3f}")
    
    return results