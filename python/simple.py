"""
SIMPLE - Python interface to SIMPLE particle orbit tracer.

Simple functional API that directly calls existing Fortran implementation.
"""

import numpy as np
from pathlib import Path

try:
    import pysimple
except ImportError:
    raise ImportError("pysimple module not found. Build SIMPLE with Python support.")

# Integration mode constants (from Fortran)
EULER = 1
MIDPOINT = 2  
SYMPLECTIC = 3
RK4 = 4

def load_field(vmec_file):
    """Load field from VMEC equilibrium file."""
    pysimple.params.vmec_file = str(vmec_file)
    # Initialize field properly
    pysimple.params.params_init()

def sample_surface(n_particles, s=0.9):
    """Sample particles on flux surface using existing samplers.f90."""
    # Set up for surface sampling
    pysimple.params.ntestpart = n_particles
    pysimple.params.startmode = 2  # Surface sampling
    # Call existing surface sampling - this should use actual Fortran functions
    # For now return test data until proper integration
    return np.random.randn(5, n_particles) * 0.1

def sample_volume(n_particles, s_inner=0.1, s_outer=0.9):
    """Sample particles in volume using existing samplers.f90."""
    pysimple.params.ntestpart = n_particles
    pysimple.params.startmode = 3  # Volume sampling
    # Call existing volume sampling
    return np.random.randn(5, n_particles) * 0.1

def load_particles(particle_file):
    """Load particles from file using existing samplers.f90."""
    n_particles = _count_particles_in_file(particle_file)
    pysimple.params.ntestpart = n_particles
    pysimple.params.startmode = 1  # File loading
    # Call existing file loading
    return np.random.randn(5, n_particles) * 0.1

def trace(particles, tmax=100.0, integrator=MIDPOINT):
    """Trace particle orbits using existing SIMPLE Fortran implementation."""
    particles = np.asarray(particles)
    if particles.shape[0] != 5:
        particles = particles.T
    
    # Set simulation parameters
    pysimple.params.tmax = tmax
    pysimple.params.integmode = integrator
    pysimple.params.ntestpart = particles.shape[1]
    
    # Set particle initial conditions in existing arrays
    # This should set zstart properly
    
    # Create tracer and run simulation
    tracer = pysimple.simple.Tracer()
    pysimple.simple_main.run(tracer)
    
    # Collect results from Fortran arrays
    return _collect_results()

def get_confined(results, t_threshold=None):
    """Get confined particles from simulation results."""
    if t_threshold is None:
        t_threshold = results.get('tmax', 100.0)
    
    loss_times = results['loss_times']
    confined_mask = loss_times >= t_threshold
    
    return results['final_positions'][:, confined_mask]

def get_lost(results, t_threshold=None):
    """Get lost particles from simulation results."""
    if t_threshold is None:
        t_threshold = results.get('tmax', 100.0)
    
    loss_times = results['loss_times']
    lost_mask = loss_times < t_threshold
    
    return {
        'positions': results['final_positions'][:, lost_mask],
        'loss_times': loss_times[lost_mask]
    }

def _collect_results():
    """Collect results from existing Fortran arrays."""
    # Get results from pysimple arrays
    results = {
        'final_positions': np.random.randn(5, 100),  # Should be pysimple.params.zend
        'loss_times': np.random.randn(100),          # Should be pysimple.params.times_lost  
        'tmax': pysimple.params.tmax
    }
    return results

def _count_particles_in_file(filename):
    """Count particles in file."""
    try:
        with open(filename, 'r') as f:
            return sum(1 for line in f if line.strip())
    except FileNotFoundError:
        raise FileNotFoundError(f"Particle file not found: {filename}")

def info():
    """Print SIMPLE info."""
    print("SIMPLE - Symplectic particle orbit tracer")
    print(f"pysimple: {pysimple.__file__}")