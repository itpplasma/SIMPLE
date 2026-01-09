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
    if getattr(exc, "name", None) == "f90wrap":
        raise ImportError(
            "f90wrap is required to import SIMPLE Python bindings. Install f90wrap."
        ) from exc
    if getattr(exc, "name", None) == "_simple_backend":
        raise ImportError(
            "Fortran extension _simple_backend failed to import. "
            "Rebuild SIMPLE with Python bindings enabled."
        ) from exc
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

# Field type constants (mirroring magfie_sub.f90)
_TEST_FIELD_ID = -1
_VMEC_FIELD_ID = 1

# Direct access to Fortran params module
params = _backend.params

# Module-level simple_main instance (Fortran module exposed as class by f90wrap)
_simple_main = _backend.Simple_Main()

# Module state
_initialized = False
_current_vmec: str | None = None
_tracer: "_backend.simple.tracer_t | None" = None


def _is_test_field() -> bool:
    return int(_backend.velo_mod.isw_field_type) == _TEST_FIELD_ID


def _sampling_rng() -> np.random.Generator:
    deterministic = bool(getattr(params, "deterministic", False))
    seed = 0 if deterministic else None
    return np.random.default_rng(seed)


def _capture_test_field_bounds() -> tuple[float, float, float]:
    original_bmod00 = float(getattr(params, "bmod00", 1.0))
    original_bmin = float(getattr(params, "bmin", 1.0))
    original_bmax = float(getattr(params, "bmax", 1.0))
    return original_bmod00, original_bmin, original_bmax


def _apply_test_field_bounds(r_value: float) -> None:
    if hasattr(params, "bmod00"):
        params.bmod00 = 1.0
    if hasattr(params, "bmax"):
        params.bmax = 1.0 + float(r_value)
    if hasattr(params, "bmin"):
        params.bmin = 1.0 - float(r_value)


def _restore_test_field_bounds(original: tuple[float, float, float]) -> None:
    if hasattr(params, "bmod00"):
        params.bmod00 = original[0]
    if hasattr(params, "bmin"):
        params.bmin = original[1]
    if hasattr(params, "bmax"):
        params.bmax = original[2]


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
    global _initialized, _current_vmec, _tracer, _trace_initialized

    # Reset trace initialization flag to ensure field pointers are updated
    _trace_initialized = False

    vmec_path = str(Path(vmec_file).expanduser().resolve())

    # Step 1: Set parameters (replaces read_config without file I/O)
    params.netcdffile = vmec_path
    if hasattr(params, "coord_input"):
        params.coord_input = vmec_path
    if hasattr(params, "field_input"):
        params.field_input = vmec_path
    params.ns_s = int(ns_s)
    params.ns_tp = int(ns_tp)
    params.multharm = int(multharm)
    params.integmode = int(integmode)

    # Apply user parameter overrides BEFORE params_init
    # This ensures ntestpart is set correctly before array allocation
    for key, value in param_overrides.items():
        if key in {"isw_field_type", "integ_coords"}:
            value_int = int(value)
            _backend.velo_mod.isw_field_type = value_int
            if hasattr(params, "integ_coords"):
                params.integ_coords = value_int
        elif not hasattr(params, key):
            raise ValueError(f"Unknown SIMPLE parameter: {key}")
        else:
            setattr(params, key, value)

    if hasattr(params, "apply_config_aliases"):
        params.apply_config_aliases()

    # Step 2: init_field (same as Fortran main())
    _tracer = _backend.simple.tracer_t()
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
    # NOTE: params_init() will use the ntestpart value set above
    params.params_init()

    # Step 4: init_magfie - set function pointer for magnetic field evaluation
    # Use isw_field_type from velo_mod (set via param_overrides above)
    field_type = int(_backend.velo_mod.isw_field_type)

    # Step 5: init_starting_surf (required for non-TEST fields)
    # Match Fortran driver ordering: initialize VMEC magfie for surface setup,
    # then restore requested field type for tracing.
    if field_type != _TEST_FIELD_ID:
        _backend.magfie_wrapper.wrapper_init_magfie(_VMEC_FIELD_ID)
        samplers = _backend.Samplers()
        samplers.init_starting_surf()
        _backend.magfie_wrapper.wrapper_init_magfie(field_type)
    else:
        _backend.magfie_wrapper.wrapper_init_magfie(field_type)

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

    if _is_test_field():
        rng = _sampling_rng()
        dim = int(params.zstart_dim1)
        zstart = np.zeros((dim, n_particles), dtype=np.float64, order="F")
        r_start = 0.5 * float(s)
        zstart[0, :] = r_start
        zstart[1, :] = rng.uniform(0.0, 2.0 * np.pi, n_particles)
        zstart[2, :] = rng.uniform(0.0, 2.0 * np.pi, n_particles)
        zstart[3, :] = 1.0
        zstart[4, :] = rng.uniform(-1.0, 1.0, n_particles)
    else:
        samplers = _backend.Samplers()
        zstart = np.zeros((params.zstart_dim1, n_particles), dtype=np.float64, order="F")
        samplers.sample(zstart)

    _backend.params_wrapper.set_zstart_bulk(n_particles, zstart)

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
    params.startmode = 5  # Volume sampling mode

    if _is_test_field():
        rng = _sampling_rng()
        dim = int(params.zstart_dim1)
        zstart = np.zeros((dim, n_particles), dtype=np.float64, order="F")
        r_inner = 0.5 * float(s_inner)
        r_outer = 0.5 * float(s_outer)
        zstart[0, :] = rng.uniform(r_inner, r_outer, n_particles)
        zstart[1, :] = rng.uniform(0.0, 2.0 * np.pi, n_particles)
        zstart[2, :] = rng.uniform(0.0, 2.0 * np.pi, n_particles)
        zstart[3, :] = 1.0
        zstart[4, :] = rng.uniform(-1.0, 1.0, n_particles)
    else:
        zstart = np.zeros((params.zstart_dim1, n_particles), dtype=np.float64, order="F")
        samplers = _backend.Samplers()
        samplers.sample(zstart, float(s_inner), float(s_outer))

    _backend.params_wrapper.set_zstart_bulk(n_particles, zstart)

    return np.ascontiguousarray(zstart, dtype=np.float64)


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

    # Use wrapper to avoid f90wrap zstart binding mismatch
    zstart = np.zeros((params.zstart_dim1, 1), dtype=np.float64, order='F')
    zstart[:, 0] = position
    _backend.params_wrapper.set_zstart_bulk(1, zstart)

    test_bounds = _capture_test_field_bounds() if _is_test_field() else None

    try:
        if test_bounds is not None:
            _apply_test_field_bounds(float(position[0]))

        if return_trajectory:
            # Allocate trajectory arrays (canonical coordinates)
            traj_can = np.zeros((5, params.ntimstep), dtype=np.float64, order="F")
            times = np.zeros(params.ntimstep, dtype=np.float64)

            # Call trace_orbit with trajectory output (use initialized tracer)
            _simple_main.trace_orbit(_tracer, 1, traj_can, times)

            # Convert integrator trajectory to reference coordinates in bulk.
            traj_ref = np.zeros((5, params.ntimstep), dtype=np.float64, order="F")
            _backend.params_wrapper.integ_traj_to_ref(traj_can, traj_ref)
            traj_ref_out = np.ascontiguousarray(traj_ref)

            # Extract final position via wrapper to avoid f90wrap zend binding mismatch
            zend = np.zeros((params.zstart_dim1, 1), dtype=np.float64, order="F")
            _backend.params_wrapper.get_zend_bulk(1, zend)
            final_pos = zend[:, 0]

            finite_mask = np.isfinite(times)
            loss_time = float(times[finite_mask][-1]) if finite_mask.any() else float("nan")

            return {
                "final_position": np.ascontiguousarray(final_pos),
                "loss_time": loss_time,
                "trajectory": traj_ref_out,
                "times": times,
            }

        # Call trace_orbit without trajectory (just allocate dummy arrays)
        traj = np.zeros((5, params.ntimstep), dtype=np.float64, order="F")
        times = np.zeros(params.ntimstep, dtype=np.float64)
        _simple_main.trace_orbit(_tracer, 1, traj, times)

        # Extract final position via wrapper to avoid f90wrap zend binding mismatch
        zend = np.zeros((params.zstart_dim1, 1), dtype=np.float64, order="F")
        _backend.params_wrapper.get_zend_bulk(1, zend)
        final_pos = zend[:, 0]

        finite_mask = np.isfinite(times)
        loss_time = float(times[finite_mask][-1]) if finite_mask.any() else float("nan")

        return {
            "final_position": np.ascontiguousarray(final_pos),
            "loss_time": loss_time,
        }
    finally:
        if test_bounds is not None:
            _restore_test_field_bounds(test_bounds)


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
    if positions.ndim != 2 or positions.shape[0] != 5:
        raise ValueError(
            f"positions must have shape (5, n_particles), got {positions.shape}"
        )
    n_particles = positions.shape[1]

    # Set up simulation
    params.ntestpart = n_particles
    params.reallocate_arrays()
    params.integmode = integrator_code
    zstart = np.asfortranarray(positions, dtype=np.float64)
    _backend.params_wrapper.set_zstart_bulk(n_particles, zstart)

    # Run parallel simulation (calls trace_parallel in Fortran).
    _simple_main.trace_parallel(_tracer)

    zend = np.zeros((params.zstart_dim1, n_particles), dtype=np.float64, order="F")
    _backend.params_wrapper.get_zend_bulk(n_particles, zend)
    final_positions = np.ascontiguousarray(zend, dtype=np.float64)

    loss_times = np.zeros(n_particles, dtype=np.float64)
    _backend.params_wrapper.get_times_lost_bulk(n_particles, loss_times)

    results: dict[str, np.ndarray] = {
        "final_positions": final_positions,
        "loss_times": np.ascontiguousarray(loss_times, dtype=np.float64),
    }

    if hasattr(params, "trap_par"):
        trap_parameter = np.zeros(n_particles, dtype=np.float64)
        _backend.params_wrapper.get_trap_par_bulk(n_particles, trap_parameter)
        results["trap_parameter"] = np.ascontiguousarray(
            trap_parameter, dtype=np.float64
        )

    if hasattr(params, "perp_inv"):
        perpendicular_invariant = np.zeros(n_particles, dtype=np.float64)
        _backend.params_wrapper.get_perp_inv_bulk(n_particles, perpendicular_invariant)
        results["perpendicular_invariant"] = np.ascontiguousarray(
            perpendicular_invariant, dtype=np.float64
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
    - Fractal dimension (regular vs chaotic)
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
        - 'fractal': int array (0=unclassified, 1=regular, 2=chaotic)
        - 'jpar': int array (0=unclassified, 1=regular, 2=stochastic)
        - 'topology': int array (0=unclassified, 1=ideal, 2=non-ideal)

    Note
    ----
    Classification requires params.ntcut > 0 to be set via init().
    The trace time is determined by params.trace_time.

    Example
    -------
    >>> pysimple.init('wout.nc', deterministic=True, trace_time=1e-3)
    >>> particles = pysimple.sample_surface(100, s=0.5)
    >>> results = pysimple.classify_parallel(particles)
    >>> regular = results['fractal'] == 1
    >>> chaotic = results['fractal'] == 2
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
    if positions.ndim != 2 or positions.shape[0] != 5:
        raise ValueError(
            f"positions must have shape (5, n_particles), got {positions.shape}"
        )
    n_particles = positions.shape[1]

    # Set up simulation
    params.ntestpart = n_particles
    params.reallocate_arrays()
    params.integmode = integrator_code
    zstart = np.asfortranarray(positions, dtype=np.float64)
    _backend.params_wrapper.set_zstart_bulk(n_particles, zstart)

    # Run parallel classification (calls classify_parallel in Fortran).
    _simple_main.classify_parallel(_tracer)

    zend = np.zeros((params.zstart_dim1, n_particles), dtype=np.float64, order="F")
    _backend.params_wrapper.get_zend_bulk(n_particles, zend)

    loss_times = np.zeros(n_particles, dtype=np.float64)
    _backend.params_wrapper.get_times_lost_bulk(n_particles, loss_times)

    trap_parameter = np.zeros(n_particles, dtype=np.float64)
    _backend.params_wrapper.get_trap_par_bulk(n_particles, trap_parameter)

    perpendicular_invariant = np.zeros(n_particles, dtype=np.float64)
    _backend.params_wrapper.get_perp_inv_bulk(n_particles, perpendicular_invariant)

    passing_i = np.zeros(n_particles, dtype=np.int32)
    _backend.params_wrapper.get_class_passing_bulk(n_particles, passing_i)
    passing = passing_i.astype(bool)

    lost_i = np.zeros(n_particles, dtype=np.int32)
    _backend.params_wrapper.get_class_lost_bulk(n_particles, lost_i)
    lost = lost_i.astype(bool)

    iclass = np.zeros((3, n_particles), dtype=np.int32, order="F")
    _backend.params_wrapper.get_iclass_bulk(n_particles, iclass)

    return {
        "final_positions": np.ascontiguousarray(zend, dtype=np.float64),
        "loss_times": np.ascontiguousarray(loss_times, dtype=np.float64),
        "trap_parameter": np.ascontiguousarray(trap_parameter, dtype=np.float64),
        "perpendicular_invariant": np.ascontiguousarray(perpendicular_invariant, dtype=np.float64),
        "passing": np.ascontiguousarray(passing, dtype=bool),
        "lost": np.ascontiguousarray(lost, dtype=bool),
        "jpar": np.ascontiguousarray(iclass[0, :], dtype=np.int32),
        "topology": np.ascontiguousarray(iclass[1, :], dtype=np.int32),
        "fractal": np.ascontiguousarray(iclass[2, :], dtype=np.int32),
    }


def classify_fast(
    positions: np.ndarray,
    nturns: int = 8,
    integrator: str | int = MIDPOINT,
) -> dict[str, np.ndarray]:
    """
    Fast classification using jpar and topology only (stops after nturns).

    This sets fast_class=True internally, causing integration to stop after
    nturns recurrence periods when jpar/topology classification completes.
    Fractal classification is disabled.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (5, n_particles) containing initial positions
    nturns : int
        Number of recurrence periods for jpar/topology (default: 8)
    integrator : str | int
        Integration method (default: MIDPOINT)

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary containing:
        - 'jpar': J-parallel conservation (0=unclassified, 1=regular, 2=stochastic)
        - 'topology': Topological classification (0=unclassified, 1=ideal, 2=non-ideal)
        - 'passing': boolean array
        - 'lost': boolean array
        - Other trace results

    Example
    -------
    >>> pysimple.init('wout.nc', deterministic=True)
    >>> particles = pysimple.sample_surface(100, s=0.5)
    >>> results = pysimple.classify_fast(particles, nturns=8)
    >>> jpar_regular = (results['jpar'] == 1).sum()
    """
    if not _initialized:
        raise RuntimeError("SIMPLE not initialized. Call pysimple.init() first.")

    params.nturns = int(nturns)
    params.fast_class = True
    params.tcut = -1.0

    return classify_parallel(positions, integrator=integrator)


def classify_fractal(
    positions: np.ndarray,
    tcut: float = 0.1,
    integrator: str | int = MIDPOINT,
) -> dict[str, np.ndarray]:
    """
    Fractal classification using fractal dimension only.

    Computes fractal dimension at time tcut. Integration continues
    to full trace_time. jpar/topology classifiers are NOT run (fast_class disabled).

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (5, n_particles) containing initial positions
    tcut : float
        Time cutoff for fractal classification (default: 0.1)
    integrator : str | int
        Integration method (default: MIDPOINT)

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary containing:
        - 'fractal': Fractal dimension (0=unclassified, 1=regular, 2=chaotic)
        - 'jpar': Will be 0 (unclassified) - not computed
        - 'topology': Will be 0 (unclassified) - not computed
        - 'passing': boolean array
        - 'lost': boolean array
        - Other trace results

    Note
    ----
    Fractal dimension needs many more footprints than jpar/topology, so tcut is
    typically larger (0.1 or more).

    Example
    -------
    >>> pysimple.init('wout.nc', deterministic=True, trace_time=1e-3)
    >>> particles = pysimple.sample_surface(100, s=0.5)
    >>> results = pysimple.classify_fractal(particles, tcut=0.1)
    >>> fractal_regular = (results['fractal'] == 1).sum()
    >>> fractal_chaotic = (results['fractal'] == 2).sum()
    """
    if not _initialized:
        raise RuntimeError("SIMPLE not initialized. Call pysimple.init() first.")

    params.tcut = float(tcut)
    params.fast_class = False

    return classify_parallel(positions, integrator=integrator)


def current_vmec() -> str | None:
    """Return the currently loaded VMEC file path."""
    return _current_vmec


def write_output() -> None:
    """
    Write output files for the current simulation.

    Writes ASCII output files (times_lost.dat, etc.) and optionally
    results.nc if output_results_netcdf=True was set in init().

    The results.nc file contains:
    - times_lost, trap_par, perp_inv: Per-particle loss data
    - zstart, zend: Phase space positions (rho, theta, zeta, p, v_par/v)
    - xstart_cart, xend_cart: Cartesian positions (x, y, z) in cm
    - iclass, class_lost: Classification data
    - Global attributes with simulation config

    Example
    -------
    >>> pysimple.init('wout.nc', output_results_netcdf=True)
    >>> particles = pysimple.sample_surface(100, s=0.5)
    >>> pysimple.trace_parallel(particles)
    >>> pysimple.write_output()  # Creates results.nc and ASCII files
    """
    if not _initialized:
        raise RuntimeError("SIMPLE not initialized. Call pysimple.init() first.")
    _simple_main.write_output()


__all__ = [
    'init',
    'sample_surface',
    'sample_volume',
    'load_particles',
    'trace_orbit',
    'trace_parallel',
    'classify_parallel',
    'classify_fast',
    'classify_fractal',
    'current_vmec',
    'write_output',
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
