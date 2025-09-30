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

# Integration mode constants (from orbit_symplectic_base.f90)
RK45 = 0
EXPL_IMPL_EULER = 1  # Default
IMPL_EXPL_EULER = 2
MIDPOINT = 3
GAUSS1 = 4
GAUSS2 = 5
GAUSS3 = 6
GAUSS4 = 7
LOBATTO3 = 15

# Default values (from params.f90)
DEFAULT_TMAX = 0.1          # trace_time=1d-1
DEFAULT_SURFACE = 0.5       # sbeg=0.5d0
DEFAULT_S_INNER = 0.1       
DEFAULT_S_OUTER = 0.9
DEFAULT_INTEGRATOR = EXPL_IMPL_EULER  # integmode default
DEFAULT_NTESTPART = 1024    # ntestpart=1024
DEFAULT_NTIMSTEP = 10000    # ntimstep=10000

def load_field(vmec_file, ns_s=5, ns_tp=5, multharm=5):
    """Load field from VMEC equilibrium file.

    This follows the initialization sequence from app/simple.f90:
    1. Set netcdf file path
    2. Initialize field (init_field)
    3. Initialize parameters (params_init)

    Args:
        vmec_file: Path to VMEC NetCDF file
        ns_s: Spline order for 3D quantities over s variable (default: 5)
        ns_tp: Spline order for 3D quantities over theta and phi (default: 5)
        multharm: Angular grid factor (default: 5)
    """
    # Create a tracer object for field initialization
    tracer = pysimple.simple.Tracer()

    # Set VMEC file path
    pysimple.params.netcdffile = str(vmec_file)

    # Initialize field with VMEC file (equivalent to init_field call in app/simple.f90:29)
    # Uses standard defaults: ns_s=5, ns_tp=5, multharm=5 (from examples/simple_full.in)
    pysimple.simple_main.init_field(
        tracer,
        str(vmec_file),
        ns_s,
        ns_tp,
        multharm,
        pysimple.params.integmode
    )

    # Initialize parameters (equivalent to params_init call in app/simple.f90:32)
    pysimple.params.params_init()

def sample_surface(n_particles, s=DEFAULT_SURFACE):
    """Sample particles on flux surface using existing samplers.f90."""
    pysimple.params.ntestpart = n_particles
    pysimple.params.startmode = 2  # Surface sampling
    pysimple.params.sbeg[0] = s  # Set surface
    pysimple.params.reallocate_arrays()

    # Call existing surface sampling through Samplers class
    samplers = pysimple.Samplers()
    samplers.sample_surface_fieldline(pysimple.params.zstart)

    # Return copy of zstart array (only the requested number of particles)
    return np.copy(pysimple.params.zstart[:, :n_particles])

def sample_volume(n_particles, s_inner=DEFAULT_S_INNER, s_outer=DEFAULT_S_OUTER):
    """Sample particles in volume using existing samplers.f90."""
    pysimple.params.ntestpart = n_particles
    pysimple.params.startmode = 3  # Volume sampling
    pysimple.params.reallocate_arrays()

    # Call existing volume sampling through Samplers class
    samplers = pysimple.Samplers()
    samplers.sample_volume_single(pysimple.params.zstart, s_inner, s_outer)

    # Return copy of zstart array (only the requested number of particles)
    return np.copy(pysimple.params.zstart[:, :n_particles])

def load_particles(particle_file):
    """Load particles from file using existing samplers.f90."""
    n_particles = _count_particles_in_file(str(particle_file))
    pysimple.params.ntestpart = n_particles
    pysimple.params.startmode = 1  # File loading
    pysimple.params.reallocate_arrays()

    # Call existing file loading through Samplers class
    samplers = pysimple.Samplers()
    samplers.sample_read(pysimple.params.zstart, str(particle_file))

    # Return copy of zstart array (only the requested number of particles)
    return np.copy(pysimple.params.zstart[:, :n_particles])

def trace(particles, tmax=DEFAULT_TMAX, integrator=DEFAULT_INTEGRATOR):
    """Trace particle orbits using existing SIMPLE Fortran implementation."""
    particles = np.asarray(particles)
    if particles.shape[0] != 5:
        particles = particles.T
    
    # Set simulation parameters
    pysimple.params.trace_time = tmax
    pysimple.params.integmode = integrator
    pysimple.params.ntestpart = particles.shape[1]
    pysimple.params.reallocate_arrays()
    
    # Set particle initial conditions
    pysimple.params.zstart[...] = particles
    
    # Create tracer and run simulation
    tracer = pysimple.simple.Tracer()
    pysimple.simple_main.run(tracer)
    
    # Collect results from Fortran arrays
    return _collect_results()

def get_confined(results, t_threshold=None):
    """Get confined particles from simulation results."""
    if t_threshold is None:
        t_threshold = results.get('tmax', DEFAULT_TMAX)
    
    loss_times = results['loss_times']
    confined_mask = loss_times >= t_threshold
    
    return results['final_positions'][:, confined_mask]

def get_lost(results, t_threshold=None):
    """Get lost particles from simulation results."""
    if t_threshold is None:
        t_threshold = results.get('tmax', DEFAULT_TMAX)
    
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
        'final_positions': np.copy(pysimple.params.zend),
        'loss_times': np.copy(pysimple.params.times_lost),
        'tmax': pysimple.params.trace_time
    }
    return results

def _count_particles_in_file(filename):
    """Count particles in file."""
    with open(filename, 'r') as f:
        count = 0
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):  # Skip comments
                count += 1
        return count

def set_parameters(**kwargs):
    """Set SIMPLE parameters directly.
    
    Common parameters:
    - ntimstep: Number of time steps (default 10000)
    - relerr: Relative error tolerance (default 1e-13)
    - facE_al: Energy factor (default 1.0)
    - n_e: Charge number (default 2)
    - n_d: Mass number (default 4)
    - deterministic: Use deterministic random seed (default False)
    - debug: Enable debug output (default False)
    """
    for key, value in kwargs.items():
        if hasattr(pysimple.params, key):
            setattr(pysimple.params, key, value)
        else:
            raise ValueError(f"Unknown parameter: {key}")

def get_parameters(*args):
    """Get SIMPLE parameters.
    
    Args:
        *args: Parameter names to get. If none given, returns dict of common parameters.
    """
    if not args:
        # Return common parameters
        common_params = ['trace_time', 'ntestpart', 'ntimstep', 'integmode', 'relerr', 
                        'facE_al', 'n_e', 'n_d', 'deterministic', 'debug']
        return {param: getattr(pysimple.params, param) for param in common_params 
                if hasattr(pysimple.params, param)}
    else:
        # Return specific parameters
        result = {}
        for param in args:
            if hasattr(pysimple.params, param):
                result[param] = getattr(pysimple.params, param)
            else:
                raise ValueError(f"Unknown parameter: {param}")
        return result if len(args) > 1 else result[args[0]]

def info():
    """Print SIMPLE info."""
    print("SIMPLE - Symplectic particle orbit tracer")
    print(f"pysimple: {pysimple.__file__}")
    print(f"Integration modes: RK45={RK45}, EXPL_IMPL_EULER={EXPL_IMPL_EULER}, MIDPOINT={MIDPOINT}, etc.")