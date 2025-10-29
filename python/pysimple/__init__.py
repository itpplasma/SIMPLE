"""
Simple module-level API for SIMPLE - mirroring Fortran global state.

Usage:
    import pysimple

    pysimple.init('wout.nc', deterministic=True, ntestpart=100)
    particles = pysimple.sample_surface(100, s=0.5)
    pysimple.trace(particles, tmax=1e-3)

    # Direct parameter access
    pysimple.params.trace_time = 2e-3
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

try:
    import simple_backend as _backend
except ImportError as exc:
    raise ImportError(
        "simple_backend module not found. Build SIMPLE with Python bindings enabled."
    ) from exc

# Integration mode constants (mirroring orbit_symplectic_base.f90)
RK45 = 0
EXPL_IMPL_EULER = 1
IMPL_EXPL_EULER = 2
MIDPOINT = 3
GAUSS1 = 4
GAUSS2 = 5
GAUSS3 = 6
GAUSS4 = 7
LOBATTO3 = 15

_INTEGRATOR_ALIASES: dict[str, int] = {
    "rk45": RK45,
    "expl_impl_euler": EXPL_IMPL_EULER,
    "explicit_implicit_euler": EXPL_IMPL_EULER,
    "impl_expl_euler": IMPL_EXPL_EULER,
    "implicit_explicit_euler": IMPL_EXPL_EULER,
    "midpoint": MIDPOINT,
    "symplectic_midpoint": MIDPOINT,
    "gauss1": GAUSS1,
    "gauss2": GAUSS2,
    "gauss3": GAUSS3,
    "gauss4": GAUSS4,
    "lobatto3": LOBATTO3,
}

# Default spline orders used by Fortran
DEFAULT_NS_S = 5
DEFAULT_NS_TP = 5
DEFAULT_MULTHARM = 5

# Direct access to Fortran params module
params = _backend.params

# Module-level simple_main instance (Fortran module exposed as class by f90wrap)
_simple_main = _backend.Simple_Main()

# Module state
_initialized = False
_current_vmec: str | None = None
_tracer: _backend.simple.Tracer | None = None


def init(
    vmec_file: str | Path,
    *,
    ns_s: int = DEFAULT_NS_S,
    ns_tp: int = DEFAULT_NS_TP,
    multharm: int = DEFAULT_MULTHARM,
    integmode: int = MIDPOINT,
    **param_overrides: Any,
) -> None:
    """
    Initialize SIMPLE with VMEC file and parameters.

    Follows the same initialization sequence as simple_main::main():
    1. Set parameters (replaces read_config, without file I/O)
    2. init_field() - loads VMEC, initializes field evaluation
    3. params_init() - computes derived parameters, resets seed if deterministic

    Parameters
    ----------
    vmec_file : str | Path
        Path to VMEC equilibrium NetCDF file
    ns_s : int
        Spline order for radial coordinate (default: 5)
    ns_tp : int
        Spline order for theta/phi coordinates (default: 5)
    multharm : int
        Multiharmonic order (default: 5)
    integmode : int
        Integration mode (default: MIDPOINT)
    **param_overrides
        Additional Fortran parameters to set (e.g., deterministic=True, ntestpart=100)

    Example
    -------
    >>> import pysimple
    >>> pysimple.init('wout.nc', deterministic=True, ntestpart=1000, trace_time=1e-3)
    """
    global _initialized, _current_vmec, _tracer

    vmec_path = str(Path(vmec_file).expanduser().resolve())

    # Step 1: Set parameters (replaces read_config without file I/O)
    params.netcdffile = vmec_path
    params.ns_s = int(ns_s)
    params.ns_tp = int(ns_tp)
    params.multharm = int(multharm)
    params.integmode = int(integmode)

    # Apply user parameter overrides
    for key, value in param_overrides.items():
        # Handle isw_field_type specially - it's in velo_mod, not params
        if key == 'isw_field_type':
            _backend.velo_mod.isw_field_type = int(value)
        elif not hasattr(params, key):
            raise ValueError(f"Unknown SIMPLE parameter: {key}")
        else:
            setattr(params, key, value)

    # Step 2: init_field (same as Fortran main())
    _tracer = _backend.simple.Tracer()
    _simple_main.init_field(
        _tracer,
        vmec_path,
        params.ns_s,
        params.ns_tp,
        params.multharm,
        params.integmode,
    )

    # Step 3: params_init (same as Fortran main())
    # This calls reset_seed_if_deterministic() internally!
    # Also calls reallocate_arrays() which allocates xstart, volstart needed by init_starting_surf
    params.params_init()

    # Step 4: init_magfie - set function pointer for magnetic field evaluation
    # Use isw_field_type from velo_mod (set via param_overrides above)
    field_type = int(_backend.velo_mod.isw_field_type)
    _backend.magfie_wrapper.wrapper_init_magfie(field_type)

    # Step 5: init_starting_surf (MUST be called before sampling particles!)
    # This integrates the magnetic field line to compute bmin, bmax
    samplers = _backend.Samplers()
    samplers.init_starting_surf()

    _initialized = True
    _current_vmec = vmec_path


def sample_surface(n_particles: int, s: float) -> np.ndarray:
    """
    Sample particles on a flux surface.

    Parameters
    ----------
    n_particles : int
        Number of particles to sample
    s : float
        Normalized toroidal flux coordinate (0 to 1)

    Returns
    -------
    np.ndarray
        Array of shape (5, n_particles) containing particle positions

    Example
    -------
    >>> particles = pysimple.sample_surface(100, s=0.5)
    """
    if not _initialized:
        raise RuntimeError("SIMPLE not initialized. Call pysimple.init() first.")

    params.ntestpart = int(n_particles)
    params.reallocate_arrays()
    params.startmode = 2

    # Use wrapper to avoid f90wrap array bug
    _backend.params_wrapper.set_sbeg(1, float(s))

    samplers = _backend.Samplers()
    # Sample directly into params.zstart (using wrapper for bulk access)
    samplers.sample(params.zstart)

    # Get results using wrapper to avoid array access bug
    zstart = np.zeros((params.zstart_dim1, n_particles), dtype=np.float64, order='F')
    _backend.params_wrapper.get_zstart_bulk(n_particles, zstart)

    return np.ascontiguousarray(zstart, dtype=np.float64)


def sample_volume(n_particles: int, s_inner: float, s_outer: float) -> np.ndarray:
    """
    Sample particles in a volume between two flux surfaces.

    Parameters
    ----------
    n_particles : int
        Number of particles to sample
    s_inner : float
        Inner flux surface (0 to 1)
    s_outer : float
        Outer flux surface (0 to 1)

    Returns
    -------
    np.ndarray
        Array of shape (5, n_particles) containing particle positions

    Example
    -------
    >>> particles = pysimple.sample_volume(1000, s_inner=0.1, s_outer=0.9)
    """
    if not _initialized:
        raise RuntimeError("SIMPLE not initialized. Call pysimple.init() first.")

    params.ntestpart = int(n_particles)
    params.reallocate_arrays()
    params.startmode = 3

    samplers = _backend.Samplers()
    samplers.sample_volume_single(params.zstart, float(s_inner), float(s_outer))

    return np.ascontiguousarray(params.zstart[:, :n_particles], dtype=np.float64)


def load_particles(particle_file: str | Path) -> np.ndarray:
    """
    Load particles from a text file.

    Parameters
    ----------
    particle_file : str | Path
        Path to particle file (e.g., start.dat)

    Returns
    -------
    np.ndarray
        Array of shape (5, n_particles) containing particle positions

    Example
    -------
    >>> particles = pysimple.load_particles('start.dat')
    """
    if not _initialized:
        raise RuntimeError("SIMPLE not initialized. Call pysimple.init() first.")

    particle_path = Path(particle_file).expanduser().resolve()

    # Read file directly in Python to avoid Fortran hardcoded 'start.dat' path
    particles_list = []
    with particle_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                values = [float(x) for x in line.split()]
                if len(values) != 5:
                    raise ValueError(f"Expected 5 values per line, got {len(values)}")
                particles_list.append(values)

    if len(particles_list) == 0:
        return np.zeros((5, 0), dtype=np.float64, order="C")

    # Convert to (5, n_particles) array
    particles = np.array(particles_list, dtype=np.float64).T

    return np.ascontiguousarray(particles, dtype=np.float64)


_trace_initialized = False

def _init_before_trace():
    """Initialize components needed before tracing (call once before first trace).

    This mimics lines 78-81 from simple_main::main():
    - init_magfie(isw_field_type) - sets magfie pointer to correct field type
    - init_counters - resets counters

    Note: init_starting_surf is called during init(), not here!
    """
    global _trace_initialized

    if _trace_initialized:
        return

    # init_magfie(isw_field_type) - set field evaluation to configured type
    _backend.magfie_wrapper.wrapper_init_magfie(_backend.velo_mod.isw_field_type)

    # init_counters - reset counters
    _simple_main.init_counters()

    _trace_initialized = True


def trace_orbit(
    position: np.ndarray,
    integrator: str | int = MIDPOINT,
    return_trajectory: bool = False,
) -> dict[str, np.ndarray]:
    """
    Trace a single particle orbit.

    Parameters
    ----------
    position : np.ndarray
        Array of shape (5,) containing initial position
    integrator : str | int
        Integration method (default: MIDPOINT)
    return_trajectory : bool
        If True, return full trajectory arrays

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary containing:
        - 'final_position': array of shape (5,)
        - 'loss_time': float
        - 'trajectory': array of shape (5, ntimstep) if return_trajectory=True
        - 'times': array of shape (ntimstep,) if return_trajectory=True

    Note
    ----
    The trace time is determined by the `trace_time` parameter passed to `init()`.
    To change trace time, call `init()` again with a different `trace_time` value.

    Example
    -------
    >>> pysimple.init('wout.nc', trace_time=1e-3)
    >>> result = pysimple.trace_orbit(position, return_trajectory=True)
    >>> trajectory = result['trajectory']
    """
    if not _initialized:
        raise RuntimeError("SIMPLE not initialized. Call pysimple.init() first.")

    # Initialize components needed before first trace
    _init_before_trace()

    # Resolve integrator
    if isinstance(integrator, str):
        key = integrator.lower()
        if key not in _INTEGRATOR_ALIASES:
            raise ValueError(f"Unknown integrator: {integrator}")
        integrator_code = _INTEGRATOR_ALIASES[key]
    else:
        integrator_code = int(integrator)

    position = np.ascontiguousarray(position, dtype=np.float64)
    if position.shape != (5,):
        raise ValueError(f"position must have shape (5,), got {position.shape}")

    # Set up simulation for single particle
    params.ntestpart = 1
    params.reallocate_arrays()
    params.integmode = integrator_code
    params.zstart[:, 0] = position

    if return_trajectory:
        # Allocate trajectory arrays (canonical coordinates)
        traj_can = np.zeros((5, params.ntimstep), dtype=np.float64, order='F')
        times = np.zeros(params.ntimstep, dtype=np.float64)

        # Call trace_orbit with trajectory output (use initialized tracer)
        _simple_main.trace_orbit(_tracer, 1, traj_can, times)

        # Convert canonical to reference coordinates (match Fortran NetCDF output)
        traj_ref = np.zeros((5, params.ntimstep), dtype=np.float64, order='C')
        xref = np.zeros(3, dtype=np.float64)

        for it in range(params.ntimstep):
            _backend.field_can_mod.can_to_ref_wrapper(traj_can[0:3, it], xref)
            traj_ref[0, it] = xref[0]  # s
            traj_ref[1, it] = xref[1]  # theta
            traj_ref[2, it] = xref[2]  # phi
            traj_ref[3, it] = traj_can[3, it]  # p_abs
            traj_ref[4, it] = traj_can[4, it]  # v_par

        return {
            'final_position': np.ascontiguousarray(params.zend[:, 0]),
            'loss_time': float(params.times_lost[0]),
            'trajectory': traj_ref,
            'times': times,
        }
    else:
        # Call trace_orbit without trajectory (just allocate dummy arrays)
        traj = np.zeros((5, params.ntimstep), dtype=np.float64, order='F')
        times = np.zeros(params.ntimstep, dtype=np.float64)
        _simple_main.trace_orbit(_tracer, 1, traj, times)

        return {
            'final_position': np.ascontiguousarray(params.zend[:, 0]),
            'loss_time': float(params.times_lost[0]),
        }


def trace_parallel(
    positions: np.ndarray,
    integrator: str | int = MIDPOINT,
) -> dict[str, np.ndarray]:
    """
    Trace multiple particle orbits in parallel (no trajectory output).

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (5, n_particles) containing initial positions
    integrator : str | int
        Integration method (default: MIDPOINT)

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary containing:
        - 'final_positions': array of shape (5, n_particles)
        - 'loss_times': array of shape (n_particles,)
        - 'trap_parameter': array of shape (n_particles,) if available
        - 'perpendicular_invariant': array of shape (n_particles,) if available

    Note
    ----
    The trace time is determined by the `trace_time` parameter passed to `init()`.
    To change trace time, call `init()` again with a different `trace_time` value.

    Example
    -------
    >>> pysimple.init('wout.nc', trace_time=1e-3)
    >>> results = pysimple.trace_parallel(particles)
    >>> lost_mask = results['loss_times'] < pysimple.params.trace_time

    Notes
    -----
    For trajectory output, use trace_orbit() for single particles or loop over particles.
    This function uses the parallel Fortran implementation for speed.
    """
    if not _initialized:
        raise RuntimeError("SIMPLE not initialized. Call pysimple.init() first.")

    # Initialize components needed before first trace
    _init_before_trace()

    # Resolve integrator
    if isinstance(integrator, str):
        key = integrator.lower()
        if key not in _INTEGRATOR_ALIASES:
            raise ValueError(f"Unknown integrator: {integrator}")
        integrator_code = _INTEGRATOR_ALIASES[key]
    else:
        integrator_code = int(integrator)

    positions = np.ascontiguousarray(positions, dtype=np.float64)
    n_particles = positions.shape[1]

    # Set up simulation
    params.ntestpart = n_particles
    params.reallocate_arrays()
    params.integmode = integrator_code
    params.zstart[:, :n_particles] = positions

    # Run parallel simulation (calls trace_parallel in Fortran, uses initialized tracer)
    _simple_main.trace_parallel(_tracer)
    # Collect results
    results = {
        'final_positions': np.ascontiguousarray(
            params.zend[:, :n_particles], dtype=np.float64
        ),
        'loss_times': np.ascontiguousarray(
            params.times_lost[:n_particles], dtype=np.float64
        ),
    }

    if hasattr(params, 'trap_par'):
        results['trap_parameter'] = np.ascontiguousarray(
            params.trap_par[:n_particles], dtype=np.float64
        )

    if hasattr(params, 'perp_inv'):
        results['perpendicular_invariant'] = np.ascontiguousarray(
            params.perp_inv[:n_particles], dtype=np.float64
        )

    return results


def classify_parallel(
    positions: np.ndarray,
    integrator: str | int = MIDPOINT,
) -> dict[str, np.ndarray]:
    """
    Trace and classify multiple particle orbits in parallel.

    Performs orbit classification including:
    - Trapped/passing determination
    - Minkowski fractal dimension (regular vs chaotic)
    - J_parallel conservation (trapped only)
    - Topological classification (trapped only)

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (5, n_particles) containing initial positions
    integrator : str | int
        Integration method (default: MIDPOINT)

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary containing:
        - 'final_positions': array of shape (5, n_particles)
        - 'loss_times': array of shape (n_particles,)
        - 'trap_parameter': array of shape (n_particles,)
        - 'perpendicular_invariant': array of shape (n_particles,)
        - 'passing': boolean array (True=passing, False=trapped)
        - 'lost': boolean array (True=lost, False=confined)
        - 'minkowski': int array (0=unclassified, 1=regular, 2=chaotic)
        - 'jpar': int array (0=unclassified, 1=regular, 2=stochastic)
        - 'topology': int array (0=unclassified, 1=ideal, 2=non-ideal)

    Note
    ----
    Classification requires params.ntcut > 0 to be set via init().
    The trace time is determined by params.trace_time.

    Example
    -------
    >>> pysimple.init('wout.nc', ntcut=3, deterministic=True, trace_time=1e-3)
    >>> particles = pysimple.sample_surface(100, s=0.5)
    >>> results = pysimple.classify_parallel(particles)
    >>> regular = results['minkowski'] == 1
    >>> chaotic = results['minkowski'] == 2
    """
    if not _initialized:
        raise RuntimeError("SIMPLE not initialized. Call pysimple.init() first.")

    # Initialize components needed before first trace
    _init_before_trace()

    # Resolve integrator
    if isinstance(integrator, str):
        key = integrator.lower()
        if key not in _INTEGRATOR_ALIASES:
            raise ValueError(f"Unknown integrator: {integrator}")
        integrator_code = _INTEGRATOR_ALIASES[key]
    else:
        integrator_code = int(integrator)

    positions = np.ascontiguousarray(positions, dtype=np.float64)
    n_particles = positions.shape[1]

    # Set up simulation
    params.ntestpart = n_particles
    params.reallocate_arrays()
    params.integmode = integrator_code
    params.zstart[:, :n_particles] = positions

    # Run parallel classification (calls classify_parallel in Fortran)
    _simple_main.classify_parallel(_tracer)

    # Collect results
    results = {
        'final_positions': np.ascontiguousarray(
            params.zend[:, :n_particles], dtype=np.float64
        ),
        'loss_times': np.ascontiguousarray(
            params.times_lost[:n_particles], dtype=np.float64
        ),
        'trap_parameter': np.ascontiguousarray(
            params.trap_par[:n_particles], dtype=np.float64
        ),
        'perpendicular_invariant': np.ascontiguousarray(
            params.perp_inv[:n_particles], dtype=np.float64
        ),
        'passing': np.ascontiguousarray(
            params.class_passing[:n_particles], dtype=bool
        ),
        'lost': np.ascontiguousarray(
            params.class_lost[:n_particles], dtype=bool
        ),
        'jpar': np.ascontiguousarray(
            params.iclass[0, :n_particles], dtype=np.int32
        ),
        'topology': np.ascontiguousarray(
            params.iclass[1, :n_particles], dtype=np.int32
        ),
        'minkowski': np.ascontiguousarray(
            params.iclass[2, :n_particles], dtype=np.int32
        ),
    }

    return results


def current_vmec() -> str | None:
    """Return the currently loaded VMEC file path."""
    return _current_vmec


__all__ = [
    'init',
    'sample_surface',
    'sample_volume',
    'load_particles',
    'trace_orbit',
    'trace_parallel',
    'classify_parallel',
    'current_vmec',
    'params',
    'RK45',
    'EXPL_IMPL_EULER',
    'IMPL_EXPL_EULER',
    'MIDPOINT',
    'GAUSS1',
    'GAUSS2',
    'GAUSS3',
    'GAUSS4',
    'LOBATTO3',
]
