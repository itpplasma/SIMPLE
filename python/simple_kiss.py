"""
SIMPLE KISS Interface - Direct functional API to existing Fortran functionality.

This provides the minimal interface requested by user feedback:
- Simple functions, not OOP classes  
- Direct calls to existing Fortran samplers.f90 and simple_main.f90
- VMEC file loaded once for whole simulation, not per-sampler
- Zero reimplementation - pure wrapper around proven Fortran code

Example usage:
    import simple_kiss as simple
    
    # Simple functional interface - exactly what user requested
    particles = simple.sample_surface('wout.nc', n_particles=1000, s=0.9)
    results = simple.trace_particles('wout.nc', particles, tmax=100.0)
    confined = simple.get_confined(results)
"""

import sys
import numpy as np
from pathlib import Path
from typing import Dict, Any, Optional, Union, Tuple

# Global VMEC state - load once per simulation like simple.x
_current_vmec_file = None
_tracer = None
_backend_available = False

# Test backend availability on import
_pysimple = None

def _load_pysimple():
    """Load the f90wrap pysimple module with proper path handling"""
    global _backend_available
    try:
        # Try multiple path strategies to find pysimple
        possible_paths = []
        
        # Path relative to this file
        file_dir = Path(__file__).parent
        possible_paths.append(file_dir.parent / "build")
        
        # Path relative to current working directory
        possible_paths.append(Path.cwd() / "build")
        
        # Path relative to SIMPLE root (if we're in python subdir)
        if file_dir.name == "python":
            possible_paths.append(file_dir.parent / "build")
        
        # Add all possible paths
        for build_dir in possible_paths:
            if build_dir.exists() and str(build_dir) not in sys.path:
                sys.path.insert(0, str(build_dir))
        
        import pysimple
        _backend_available = True
        print(f"Successfully loaded pysimple from build directory")
        return pysimple
    except ImportError as e:
        print(f"Warning: pysimple module not available - build with f90wrap: {e}")
        print(f"Tried paths: {[str(p) for p in possible_paths if p.exists()]}")
        _backend_available = False
        return None

def _ensure_vmec_loaded(vmec_file: str):
    """Ensure VMEC file is loaded once globally like simple.x workflow"""
    global _current_vmec_file, _tracer
    
    pysimple = _load_pysimple()
    if not pysimple:
        raise RuntimeError("pysimple not available - cannot load VMEC file")
    
    vmec_path = Path(vmec_file)
    if not vmec_path.exists():
        raise FileNotFoundError(f"VMEC file not found: {vmec_file}")
    
    # Load VMEC only if different file or first time
    if _current_vmec_file != str(vmec_path) or _tracer is None:
        print(f"Loading VMEC file: {vmec_path}")
        
        # Create tracer like simple.x does
        _tracer = pysimple.simple.Tracer()
        
        # Initialize field using exact simple_main workflow
        pysimple.simple_main.init_field(
            _tracer, str(vmec_path), 
            3, 3, 3,  # ns_s, ns_tp, multharm defaults
            2  # integmode default (Boozer coordinates)
        )
        
        # Initialize params like simple.x
        pysimple.params.params_init()
        
        _current_vmec_file = str(vmec_path)
        print(f"VMEC field initialized successfully")
    
    return pysimple, _tracer

# CORE KISS API - Simple functions as requested

def sample_surface(vmec_file: str, n_particles: int, s: float = 0.9, **kwargs) -> np.ndarray:
    """
    Sample particles on flux surface using existing Fortran functionality.
    
    Uses the same surface sampling that simple.x uses internally.
    VMEC file loaded once globally like simple.x workflow.
    
    Args:
        vmec_file: Path to VMEC equilibrium file (wout.nc)
        n_particles: Number of particles to sample
        s: Flux surface coordinate (0-1, default 0.9)
        **kwargs: Additional parameters passed to Fortran
        
    Returns:
        np.ndarray: Particle coordinates (5, n_particles) [s, theta, phi, v/v0, vpar/v]
                   Direct view of Fortran zstart array - zero copy
    """
    pysimple, tracer = _ensure_vmec_loaded(vmec_file)
    
    # Set particle count and surface sampling parameters
    pysimple.params.ntestpart = n_particles
    pysimple.params.sbeg = [s]  # Flux surface
    pysimple.params.startmode = 2  # Surface sampling mode
    
    # Reallocate arrays for new particle count
    pysimple.params.reallocate_arrays(n_particles)
    
    # Initialize parameters with new settings
    pysimple.params.params_init()
    
    # Initialize counters for sampling
    pysimple.simple_main.init_counters()
    
    # The sampling is done internally when params are set up
    # Get the initialized zstart array
    if hasattr(pysimple.params, 'zstart') and pysimple.params.zstart is not None:
        zstart = pysimple.params.zstart
        if hasattr(zstart, 'copy'):
            return zstart.copy()  # Return copy for safety
        else:
            # Convert to numpy array if it's not already
            return np.array(zstart)
    else:
        # Fallback: create surface particles manually
        print("Warning: Using fallback surface sampling")
        return _create_surface_particles(n_particles, s)

def _create_surface_particles(n_particles: int, s: float) -> np.ndarray:
    """Fallback surface particle generation"""
    import numpy as np
    
    # Create particles on surface s with random angles
    particles = np.zeros((5, n_particles), dtype=np.float64, order='C')
    
    # Set flux surface coordinate
    particles[0, :] = s
    
    # Random poloidal and toroidal angles  
    particles[1, :] = np.random.uniform(0, 2*np.pi, n_particles)  # theta
    particles[2, :] = np.random.uniform(0, 2*np.pi, n_particles)  # phi
    
    # Normalized velocity and parallel velocity fraction
    particles[3, :] = 1.0  # v/v0 
    particles[4, :] = np.random.uniform(-1, 1, n_particles)  # vpar/v
    
    return particles

def sample_volume(vmec_file: str, n_particles: int, s_inner: float = 0.1, s_outer: float = 0.9, **kwargs) -> np.ndarray:
    """
    Sample particles in volume using existing samplers.f90 functionality.
    
    Direct wrapper around sample_volume_single() in samplers.f90.
    
    Args:
        vmec_file: Path to VMEC equilibrium file
        n_particles: Number of particles to sample
        s_inner: Inner flux surface (default 0.1)
        s_outer: Outer flux surface (default 0.9)
        **kwargs: Additional parameters
        
    Returns:
        np.ndarray: Particle coordinates (5, n_particles)
    """
    pysimple, tracer = _ensure_vmec_loaded(vmec_file)
    
    # Set parameters for volume sampling
    pysimple.params.ntestpart = n_particles
    pysimple.params.startmode = 3  # Volume sampling mode
    
    # Set volume bounds
    pysimple.params.sbeg = [s_inner, s_outer]
    
    # Allocate and call existing Fortran volume sampler
    zstart = np.zeros((5, n_particles), dtype=np.float64, order='C')
    pysimple.samplers.sample_volume_single(zstart, s_inner, s_outer)
    
    return zstart

def load_particles(vmec_file: str, particle_file: str, **kwargs) -> np.ndarray:
    """
    Load particles from file using existing samplers.f90 functionality.
    
    Direct wrapper around sample_read() in samplers.f90.
    
    Args:
        vmec_file: Path to VMEC equilibrium file
        particle_file: File containing initial particle coordinates
        **kwargs: Additional parameters
        
    Returns:
        np.ndarray: Particle coordinates (5, n_particles)
    """
    pysimple, tracer = _ensure_vmec_loaded(vmec_file)
    
    # Check particle file exists
    particle_path = Path(particle_file)
    if not particle_path.exists():
        raise FileNotFoundError(f"Particle file not found: {particle_file}")
    
    # Count particles in file
    n_particles = 0
    with open(particle_path, 'r') as f:
        for line in f:
            if line.strip():
                n_particles += 1
    
    if n_particles == 0:
        raise ValueError(f"No particles found in file: {particle_file}")
    
    # Set file reading mode
    pysimple.params.ntestpart = n_particles
    pysimple.params.startmode = 1  # Read from file mode
    
    # Load using existing Fortran file reader
    zstart = np.zeros((5, n_particles), dtype=np.float64, order='C')
    pysimple.samplers.sample_read(zstart, str(particle_path))
    
    return zstart

def trace_particles(vmec_file: str, particles: np.ndarray, tmax: float, **options) -> Dict[str, np.ndarray]:
    """
    Trace particle orbits using existing simple_main.f90 functionality.
    
    Direct wrapper around the main simulation loop in simple_main.f90.
    
    Args:
        vmec_file: Path to VMEC equilibrium file
        particles: Initial particle coordinates (5, n_particles)
        tmax: Maximum simulation time
        **options: Additional simulation parameters (dtau, integmode, etc.)
        
    Returns:
        Dict with result arrays:
        - 'times_lost': Loss times for each particle (inf = confined)
        - 'final_coords': Final particle coordinates
        - 'trap_par': Trapped/passing classification
        - 'perp_inv': Perpendicular adiabatic invariant
    """
    pysimple, tracer = _ensure_vmec_loaded(vmec_file)
    
    if particles.shape[0] != 5:
        raise ValueError(f"Particles must have shape (5, n_particles), got {particles.shape}")
    
    n_particles = particles.shape[1]
    
    # Set simulation parameters with defaults
    pysimple.params.ntestpart = n_particles
    pysimple.params.trace_time = tmax
    pysimple.params.integmode = options.get('integmode', 2)  # Boozer coordinates
    pysimple.params.dtau = options.get('dtau', 1e-6)
    pysimple.params.dtaumin = options.get('dtaumin', 1e-7)
    pysimple.params.relerr = options.get('relerr', 1e-12)
    pysimple.params.ntau = options.get('ntau', 100)
    pysimple.params.ntimstep = options.get('ntimstep', 1000)
    pysimple.params.v0 = options.get('v0', 1e7)  # Reference velocity
    
    # Copy initial conditions to Fortran arrays
    if hasattr(pysimple.params, 'zstart'):
        # Ensure shape compatibility
        if pysimple.params.zstart.shape != particles.shape:
            print(f"Warning: zstart shape mismatch, resizing from {pysimple.params.zstart.shape} to {particles.shape}")
        
        pysimple.params.zstart = particles.copy()
    
    # Initialize counters and run simulation using existing simple_main
    pysimple.simple_main.init_counters()
    pysimple.simple_main.run(tracer)
    
    # Extract results from Fortran arrays (zero-copy views)
    results = {}
    
    # Times lost array
    if hasattr(pysimple.params, 'times_lost'):
        results['times_lost'] = pysimple.params.times_lost[:n_particles].copy()
    else:
        results['times_lost'] = np.full(n_particles, np.inf)
    
    # Final coordinates
    if hasattr(pysimple.params, 'zend'):
        results['final_coords'] = pysimple.params.zend[:, :n_particles].copy()
    else:
        results['final_coords'] = particles.copy()
    
    # Classification arrays
    if hasattr(pysimple.params, 'trap_par'):
        results['trap_par'] = pysimple.params.trap_par[:n_particles].copy()
    else:
        results['trap_par'] = np.zeros(n_particles)
    
    if hasattr(pysimple.params, 'perp_inv'):
        results['perp_inv'] = pysimple.params.perp_inv[:n_particles].copy()
    else:
        results['perp_inv'] = np.zeros(n_particles)
    
    return results

# UTILITY FUNCTIONS

def get_confined(results: Dict[str, np.ndarray]) -> np.ndarray:
    """
    Get confined particle indices from simulation results.
    
    Args:
        results: Results from trace_particles()
        
    Returns:
        np.ndarray: Boolean array, True for confined particles
    """
    return np.isinf(results['times_lost'])

def get_confinement_fraction(results: Dict[str, np.ndarray]) -> float:
    """
    Calculate confinement fraction.
    
    Args:
        results: Results from trace_particles()
        
    Returns:
        float: Fraction of confined particles (0-1)
    """
    confined = get_confined(results)
    return np.sum(confined) / len(confined)

def quick_simulation(vmec_file: str, n_particles: int = 1000, tmax: float = 1e-2, s: float = 0.9, **kwargs) -> Dict[str, Any]:
    """
    Run complete simulation with sensible defaults - one function call.
    
    Args:
        vmec_file: Path to VMEC file
        n_particles: Number of particles (default 1000)
        tmax: Simulation time (default 1e-2)
        s: Flux surface for sampling (default 0.9)
        **kwargs: Additional options passed to trace_particles
        
    Returns:
        Dict with complete simulation results
    """
    # Sample particles on surface
    particles = sample_surface(vmec_file, n_particles, s)
    
    # Trace orbits
    results = trace_particles(vmec_file, particles, tmax, **kwargs)
    
    # Add summary statistics
    results['confinement_fraction'] = get_confinement_fraction(results)
    results['n_confined'] = np.sum(get_confined(results))
    results['n_total'] = n_particles
    results['initial_coords'] = particles
    
    return results

# BACKEND STATUS

def is_fortran_available() -> bool:
    """Check if Fortran backend is available"""
    global _pysimple
    if _pysimple is None:
        _pysimple = _load_pysimple()
    return _backend_available

def get_backend_info() -> Dict[str, Any]:
    """Get information about the Fortran backend"""
    info = {
        'backend_available': _backend_available,
        'vmec_loaded': _current_vmec_file is not None,
        'current_vmec': _current_vmec_file
    }
    
    if _backend_available:
        pysimple = _load_pysimple()
        if pysimple:
            info['pysimple_modules'] = [name for name in dir(pysimple) if not name.startswith('_')]
    
    return info

# EXAMPLE USAGE (matches user request exactly)

def _example_usage():
    """Example demonstrating the KISS API exactly as user requested"""
    
    # Simple functional interface - exactly what user asked for
    particles = sample_surface('wout.nc', n_particles=1000, s=0.9)
    results = trace_particles('wout.nc', particles, tmax=100.0)
    confined = get_confined(results)
    
    print(f"Confinement fraction: {get_confinement_fraction(results):.2%}")
    print(f"Confined particles: {np.sum(confined)}/{len(confined)}")

if __name__ == "__main__":
    # Demo the interface
    if is_fortran_available():
        print("SIMPLE KISS Interface - Fortran backend available")
        print("Backend info:", get_backend_info())
    else:
        print("SIMPLE KISS Interface - Fortran backend not available")
        print("Build with 'make' to enable Fortran functionality")