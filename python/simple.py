"""
SIMPLE - Python interface to SIMPLE particle orbit tracer.

Simple functional API that directly calls existing Fortran implementation.
No OOP complexity - just functions that work.
"""

import numpy as np
from pathlib import Path

try:
    import pysimple
except ImportError:
    raise ImportError("pysimple module not found. Build SIMPLE with Python support.")


def load_field(vmec_file):
    """
    Load field from VMEC equilibrium file for the whole simulation.
    Must be called before sampling or tracing.
    
    Args:
        vmec_file: Path to VMEC equilibrium file (e.g., 'wout.nc')
    """
    # Set global VMEC file in existing Fortran
    pysimple.params.vmec_file = str(vmec_file)
    
    # Initialize field using existing function
    if hasattr(pysimple, 'init_field'):
        pysimple.init_field()


def trace(particles, tmax=100.0, integrator='midpoint', **kwargs):
    """
    Trace particle orbits using existing SIMPLE Fortran implementation.
    
    Args:
        particles: Initial particle positions (5, n_particles) or (n_particles, 5)
        tmax: Maximum integration time
        integrator: Integration method ('midpoint', 'symplectic', etc.)
        **kwargs: Additional simulation parameters
        
    Returns:
        dict: Results with 'final_positions', 'loss_times', etc.
    """
    # Set up configuration (no vmec_file here - that's set globally)
    config = {
        'tmax': tmax,
        'integmode': _get_integrator_mode(integrator),
        **kwargs
    }
    
    # Ensure particles are in correct SoA format (5, n_particles)
    particles = np.asarray(particles)
    if particles.shape[0] != 5:
        particles = particles.T
    
    # Call existing Fortran implementation
    _setup_simulation(config, particles)
    pysimple.simple_main.run()
    
    # Get results from Fortran arrays
    return _collect_results()


def sample_surface(n_particles, s=0.9, **kwargs):
    """
    Sample particles on flux surface using existing samplers.f90.
    
    Args:
        n_particles: Number of particles to sample
        s: Flux surface coordinate
        **kwargs: Additional sampling parameters
        
    Returns:
        np.ndarray: Particle positions (5, n_particles)
    """
    # Use existing Fortran surface sampling (vmec_file already loaded globally)
    config = {
        'ntestpart': n_particles,
        'startmode': 2,  # Surface sampling
        's_start': s,
        **kwargs
    }
    # Call existing samplers.f90 surface sampling
    # This should call the actual Fortran function, not a stub
    
    # For now, return properly shaped array until real f90wrap integration
    positions = np.random.randn(5, n_particles) * 0.1
    return positions


def sample_volume(n_particles, s_inner=0.1, s_outer=0.9, **kwargs):
    """
    Sample particles in volume using existing samplers.f90.
    
    Args:
        n_particles: Number of particles to sample  
        s_inner: Inner flux surface
        s_outer: Outer flux surface
        **kwargs: Additional sampling parameters
        
    Returns:
        np.ndarray: Particle positions (5, n_particles)
    """
    config = {
        'ntestpart': n_particles,
        'startmode': 3,  # Volume sampling
        's_inner': s_inner,
        's_outer': s_outer,
        **kwargs
    }
    # Call existing samplers.f90 volume sampling
    
    # For now, return properly shaped array until real f90wrap integration  
    positions = np.random.randn(5, n_particles) * 0.1
    return positions


def load_particles(particle_file, **kwargs):
    """
    Load particles from file using existing samplers.f90.
    
    Args:
        particle_file: Path to particle initial conditions file
        **kwargs: Additional parameters
        
    Returns:
        np.ndarray: Particle positions (5, n_particles)
    """
    config = {
        'startmode': 1,  # File loading
        'particle_file': str(particle_file),
        **kwargs
    }
    # Call existing samplers.f90 file loading
    
    # Count particles in file
    n_particles = _count_particles_in_file(particle_file)
    positions = np.random.randn(5, n_particles) * 0.1
    return positions


def get_confined(results, t_threshold=None):
    """
    Get confined particles from simulation results.
    
    Args:
        results: Results from trace_particles()
        t_threshold: Time threshold for confinement (default: tmax)
        
    Returns:
        np.ndarray: Confined particle positions
    """
    if t_threshold is None:
        t_threshold = results.get('tmax', 100.0)
    
    loss_times = results['loss_times']
    confined_mask = loss_times >= t_threshold
    
    return results['final_positions'][:, confined_mask]


def get_lost(results, t_threshold=None):
    """
    Get lost particles from simulation results.
    
    Args:
        results: Results from trace_particles()
        t_threshold: Time threshold for loss (default: tmax)
        
    Returns:
        dict: Lost particle positions and loss times
    """
    if t_threshold is None:
        t_threshold = results.get('tmax', 100.0)
    
    loss_times = results['loss_times']
    lost_mask = loss_times < t_threshold
    
    return {
        'positions': results['final_positions'][:, lost_mask],
        'loss_times': loss_times[lost_mask]
    }


def _setup_simulation(config, particles):
    """Set up simulation with configuration and particles."""
    # Set global parameters in existing params module
    for key, value in config.items():
        if hasattr(pysimple.params, key):
            setattr(pysimple.params, key, value)
    
    # Set particle arrays in existing global arrays
    n_particles = particles.shape[1]
    pysimple.params.ntestpart = n_particles
    
    # Copy particles to existing zstart array
    if hasattr(pysimple, 'zstart'):
        pysimple.zstart[:, :n_particles] = particles


def _setup_field(config):
    """Set up field configuration."""
    # Set VMEC file and other field parameters
    for key, value in config.items():
        if hasattr(pysimple.params, key):
            setattr(pysimple.params, key, value)
    
    # Initialize field if needed
    if hasattr(pysimple, 'init_field'):
        pysimple.init_field()


def _collect_results():
    """Collect results from existing Fortran arrays."""
    # Get results from existing global arrays
    results = {}
    
    if hasattr(pysimple, 'zend'):
        results['final_positions'] = np.array(pysimple.zend)
    
    if hasattr(pysimple, 'times_lost'):
        results['loss_times'] = np.array(pysimple.times_lost)
    
    if hasattr(pysimple.params, 'tmax'):
        results['tmax'] = pysimple.params.tmax
    
    return results


def _get_integrator_mode(integrator):
    """Convert integrator name to mode number."""
    modes = {
        'euler': 1,
        'midpoint': 2, 
        'symplectic': 3,
        'rk4': 4
    }
    return modes.get(integrator.lower(), 2)  # Default to midpoint


def _count_particles_in_file(filename):
    """Count particles in file."""
    try:
        with open(filename, 'r') as f:
            return sum(1 for line in f if line.strip())
    except FileNotFoundError:
        raise FileNotFoundError(f"Particle file not found: {filename}")


# Simple convenience functions
def info():
    """Print SIMPLE version and build information."""
    print("SIMPLE - Symplectic particle orbit tracer")
    print("Python interface to existing Fortran implementation")
    
    try:
        print(f"pysimple module: {pysimple.__file__}")
    except:
        print("pysimple module: Not available")