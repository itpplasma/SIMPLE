"""Python interface to SIMPLE.

pysimple is a thin layer over the Fortran module ``simple_api``. A run is
configured by the ``&config`` namelist of ``simple.in``; every keyword argument
of :func:`init` is a ``&config`` name (see ``examples/simple_full.in``)::

    import pysimple

    pysimple.init(config="simple.in")             # exactly as simple.x
    pysimple.init("wout.nc", trace_time=1e-3)     # equilibrium plus settings
    pysimple.init("wout.nc", config="simple.in")  # file with overrides

    start = pysimple.sample()                     # as configured by startmode
    result = pysimple.trace_parallel(start)

States are ``(s, theta, phi, v/v0, v_par/v)`` in reference coordinates and
have shape ``(5,)`` for one particle or ``(5, n)`` for ``n`` particles.
Angles are radians. ``pysimple.params`` exposes the Fortran ``params`` module
for reading derived values such as ``params.v0``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

try:
    import simple_backend as _backend
except ImportError as exc:
    raise ImportError(
        "SIMPLE Python bindings not found. Build SIMPLE with "
        "ENABLE_PYTHON_INTERFACE=ON or install it with pip."
    ) from exc

_api = _backend.simple_api
params = _backend.params

_integrators = _backend.orbit_symplectic_base
RK45 = _integrators.rk45
EXPL_IMPL_EULER = _integrators.expl_impl_euler
IMPL_EXPL_EULER = _integrators.impl_expl_euler
MIDPOINT = _integrators.midpoint
GAUSS1 = _integrators.gauss1
GAUSS2 = _integrators.gauss2
GAUSS3 = _integrators.gauss3
GAUSS4 = _integrators.gauss4
LOBATTO3 = _integrators.lobatto3

_INTEGRATORS = {
    "rk45": RK45,
    "expl_impl_euler": EXPL_IMPL_EULER,
    "impl_expl_euler": IMPL_EXPL_EULER,
    "midpoint": MIDPOINT,
    "gauss1": GAUSS1,
    "gauss2": GAUSS2,
    "gauss3": GAUSS3,
    "gauss4": GAUSS4,
    "lobatto3": LOBATTO3,
}

_frequencies = _backend.orbit_frequencies
ORBIT_CLASS_UNKNOWN = _frequencies.orbit_class_unknown
ORBIT_CLASS_TRAPPED = _frequencies.orbit_class_trapped
ORBIT_CLASS_PASSING = _frequencies.orbit_class_passing
FREQ_SUCCESS = _frequencies.freq_success
FREQ_INVALID_INPUT = _frequencies.freq_invalid_input
FREQ_ORBIT_LOST = _frequencies.freq_orbit_lost
FREQ_INTEGRATOR_ERROR = _frequencies.freq_integrator_error
FREQ_MAX_STEPS = _frequencies.freq_max_steps

_invariants = _backend.orbit_invariants
INVARIANT_SUCCESS = _invariants.invariant_success
INVARIANT_INVALID_INPUT = _invariants.invariant_invalid_input
INVARIANT_NOT_AXISYMMETRIC = _invariants.invariant_not_axisymmetric
INVARIANT_NO_INTERSECTION = _invariants.invariant_no_intersection
INVARIANT_ROOT_FAILURE = _invariants.invariant_root_failure
INVARIANT_CAPACITY = _invariants.invariant_capacity
SECTION_EXTREMUM_UNKNOWN = _invariants.section_extremum_unknown
SECTION_B_MINIMUM = _invariants.section_b_minimum
SECTION_B_MAXIMUM = _invariants.section_b_maximum

CUT_ANY = -1
CUT_TIP = 0
CUT_TOROIDAL = 1
_CUTS = {"any": CUT_ANY, "tip": CUT_TIP, "toroidal": CUT_TOROIDAL}

_ORBIT_CLASS_NAMES = {
    ORBIT_CLASS_UNKNOWN: "unknown",
    ORBIT_CLASS_TRAPPED: "trapped",
    ORBIT_CLASS_PASSING: "passing",
}
_SECTION_NAMES = {
    SECTION_EXTREMUM_UNKNOWN: "unknown",
    SECTION_B_MINIMUM: "B_min",
    SECTION_B_MAXIMUM: "B_max",
}
_CUT_NAMES = {CUT_TIP: "tip", CUT_TOROIDAL: "toroidal"}

_initialized = False

# Record the compiled defaults before anything changes a &config value.
_api.simple_reset_config()


def _check(status_message: tuple[int, str], what: str) -> None:
    status, message = status_message
    if status != 0:
        raise ValueError(f"{what}: {message.strip() or f'status {status}'}")


def _namelist_value(value: Any) -> str:
    if isinstance(value, (bool, np.bool_)):
        return ".true." if value else ".false."
    if isinstance(value, (str, Path)):
        return "'" + str(value).replace("'", "''") + "'"
    if np.ndim(value) > 0:
        return ", ".join(_namelist_value(item) for item in np.ravel(value))
    return repr(value.item() if isinstance(value, np.generic) else value)


def _set_config(**settings: Any) -> None:
    for name, value in settings.items():
        text = f"&config {name} = {_namelist_value(value)} /"
        _check(_api.simple_set_config(text), f"&config {name}")


def _require_init() -> None:
    if not _initialized:
        raise RuntimeError("SIMPLE not initialized. Call pysimple.init() first.")


def _integrator_code(integrator: str | int) -> int:
    if isinstance(integrator, str):
        try:
            return _INTEGRATORS[integrator.lower()]
        except KeyError:
            raise ValueError(f"Unknown integrator: {integrator}") from None
    return int(integrator)


def _use_integrator(integrator: str | int | None) -> None:
    if integrator is not None:
        params.integmode = _integrator_code(integrator)


def _state(position: np.ndarray) -> np.ndarray:
    state = np.asarray(position, dtype=np.float64)
    if state.shape != (5,):
        raise ValueError(f"position must have shape (5,), got {state.shape}")
    return state


def _states(positions: np.ndarray) -> np.ndarray:
    states = np.asfortranarray(positions, dtype=np.float64)
    if states.ndim != 2 or states.shape[0] != 5:
        raise ValueError(
            f"positions must have shape (5, n_particles), got {states.shape}"
        )
    return states


def init(
    field: str | Path | None = None,
    *,
    config: str | Path | None = None,
    plasma_profiles: np.ndarray | None = None,
    **settings: Any,
) -> None:
    """Configure SIMPLE and load the field.

    All ``&config`` values start from their defaults. ``config`` reads a
    ``simple.in`` file (its ``&config`` and ``&profiles`` groups), ``field``
    sets ``netcdffile``, ``field_input`` and ``coord_input`` to one
    equilibrium file, and ``settings`` set further ``&config`` values. Later
    sources override earlier ones in this order. A Boozer chartmap selects
    ``isw_field_type=2`` unless the field type is set.

    Raises ``ValueError`` for an unreadable file or an unknown ``&config``
    name.

    ``plasma_profiles`` is an optional ``(6, n)`` table with collision radius,
    electron/ion-1/ion-2 temperatures [eV] and ion-1/ion-2 densities [m^-3].
    Radius must increase from 0 to 1: toroidal flux ``s`` for VMEC/Boozer.
    The table is linearly interpolated, without a polynomial fit.
    """
    global _initialized

    _initialized = False
    _api.simple_reset_config()
    if config is not None:
        path = Path(config).expanduser()
        _check(_api.simple_read_config(str(path)), str(path))
    if field is not None:
        path = str(Path(field).expanduser().resolve())
        settings = {
            "netcdffile": path,
            "field_input": path,
            "coord_input": path,
            **settings,
        }
    _set_config(**settings)
    if plasma_profiles is not None:
        table = np.asfortranarray(plasma_profiles, dtype=np.float64)
        if table.ndim != 2 or table.shape[0] != 6:
            raise ValueError("plasma_profiles must have shape (6, n)")
        if _api.simple_set_profile_table(table.shape[1], table) != 0:
            raise ValueError(
                "Invalid plasma_profiles: require finite nonnegative "
                "values and strictly increasing radius from 0 to 1"
            )
    _api.simple_init()
    _initialized = True


def _coordinate_array(points: np.ndarray) -> tuple[np.ndarray, bool]:
    values = np.asarray(points, dtype=float)
    single = values.shape == (3,)
    if single:
        values = values[:, None]
    if values.ndim != 2 or values.shape[0] != 3:
        raise ValueError("coordinates must have shape (3,) or (3, n)")
    return np.asfortranarray(values), single


def cylindrical_to_reference(points: np.ndarray) -> np.ndarray:
    """Convert ``(R [m], phi [rad], Z [m])`` to loaded reference coordinates.

    Raise ValueError if an inverse fails or a physical VMEC/chartmap
    reference is unavailable (normalized TEST and SPECTRE fields).
    """
    _require_init()
    values, single = _coordinate_array(points)
    result = np.empty_like(values, order="F")
    status = np.empty(values.shape[1], dtype=np.int32)
    _api.simple_cylindrical_to_reference(values.shape[1], values, result, status)
    if np.any(status != 0):
        raise ValueError(
            f"Cylindrical inverse failed for markers "
            f"{np.flatnonzero(status != 0).tolist()}"
        )
    return result[:, 0].copy() if single else np.ascontiguousarray(result)


def reference_to_cylindrical(points: np.ndarray) -> np.ndarray:
    """Convert loaded reference coordinates to ``(R [m], phi [rad], Z [m])``."""
    _require_init()
    values, single = _coordinate_array(points)
    result = np.empty_like(values, order="F")
    status = np.empty(values.shape[1], dtype=np.int32)
    _api.simple_reference_to_cylindrical(values.shape[1], values, result, status)
    if np.any(status != 0):
        raise ValueError("Physical reference coordinates are unavailable or invalid")
    return result[:, 0].copy() if single else np.ascontiguousarray(result)


def sample(n_particles: int | None = None) -> np.ndarray:
    """Sample start states as configured by ``startmode``, ``sbeg`` and
    ``num_surf``. ``n_particles`` defaults to ``ntestpart``."""
    _require_init()
    n = int(params.ntestpart if n_particles is None else n_particles)
    z = np.zeros((5, n), order="F")
    _api.simple_sample(n, z)
    return np.ascontiguousarray(z)


def sample_surface(n_particles: int, s: float) -> np.ndarray:
    """Sample start states on flux surface ``s`` along the starting field line."""
    _require_init()
    z = np.zeros((5, n_particles), order="F")
    _api.simple_sample_surface(n_particles, float(s), z)
    return np.ascontiguousarray(z)


def sample_volume(n_particles: int, s_inner: float, s_outer: float) -> np.ndarray:
    """Sample start states uniformly in the volume ``s_inner < s < s_outer``."""
    _require_init()
    z = np.zeros((5, n_particles), order="F")
    _api.simple_sample_volume(n_particles, float(s_inner), float(s_outer), z)
    return np.ascontiguousarray(z)


def load_particles(particle_file: str | Path) -> np.ndarray:
    """Read start states from a ``start.dat``-style text file."""
    return np.loadtxt(Path(particle_file).expanduser(), ndmin=2).T.copy()


def _trace_results(n: int) -> dict[str, np.ndarray]:
    z_end = np.zeros((5, n), order="F")
    loss_times = np.zeros(n)
    trap_parameter = np.zeros(n)
    perpendicular_invariant = np.zeros(n)
    _api.simple_get_results(n, z_end, loss_times, trap_parameter,
                            perpendicular_invariant)
    return {
        "final_positions": np.ascontiguousarray(z_end),
        "loss_times": loss_times,
        "trap_parameter": trap_parameter,
        "perpendicular_invariant": perpendicular_invariant,
    }


def trace_parallel(
    positions: np.ndarray, integrator: str | int | None = None
) -> dict[str, np.ndarray]:
    """Trace particles to ``trace_time`` or loss with OpenMP, as simple.x.

    ``integrator`` overrides ``integmode``. ``loss_times`` is ``trace_time``
    for confined, ``-1`` for skipped and ``NaN`` for numerically failed
    particles.
    """
    _require_init()
    states = _states(positions)
    _use_integrator(integrator)
    _api.simple_set_start(states.shape[1], states)
    _api.simple_trace()
    return _trace_results(states.shape[1])


def classify_parallel(
    positions: np.ndarray, integrator: str | int | None = None
) -> dict[str, np.ndarray]:
    """Trace particles with the classifiers selected by ``tcut``,
    ``fast_class`` and ``nturns``.

    Adds ``passing``, ``lost`` and the ``jpar``, ``topology`` and ``fractal``
    classes (0 unclassified, 1 regular, 2 chaotic) to the trace results.
    """
    _require_init()
    states = _states(positions)
    n = states.shape[1]
    _use_integrator(integrator)
    _api.simple_set_start(n, states)
    _api.simple_classify()
    results = _trace_results(n)
    passing = np.zeros(n, dtype=np.int32)
    lost = np.zeros(n, dtype=np.int32)
    classes = np.zeros((3, n), dtype=np.int32, order="F")
    _api.simple_get_classification(n, passing, lost, classes)
    results.update(
        passing=passing.astype(bool),
        lost=lost.astype(bool),
        jpar=classes[0].copy(),
        topology=classes[1].copy(),
        fractal=classes[2].copy(),
    )
    return results


def classify_fast(
    positions: np.ndarray, nturns: int = 8, integrator: str | int | None = None
) -> dict[str, np.ndarray]:
    """Classify by J_parallel and topology, stopping after ``nturns``."""
    _require_init()
    _set_config(nturns=int(nturns), fast_class=True, tcut=-1.0)
    return classify_parallel(positions, integrator)


def classify_fractal(
    positions: np.ndarray, tcut: float = 0.1, integrator: str | int | None = None
) -> dict[str, np.ndarray]:
    """Classify by fractal dimension at time ``tcut``."""
    _require_init()
    _set_config(tcut=float(tcut), fast_class=False)
    return classify_parallel(positions, integrator)


def write_output() -> None:
    """Write the simple.x output files for the last trace."""
    _require_init()
    _api.simple_write_output()


def trace_orbit(
    position: np.ndarray,
    integrator: str | int | None = None,
    return_trajectory: bool = False,
) -> dict[str, np.ndarray | float]:
    """Trace one orbit; optionally return its ``ntimstep`` samples."""
    _require_init()
    state = _state(position)
    _use_integrator(integrator)
    nt = int(params.ntimstep)
    trajectory = np.zeros((5, nt), order="F")
    times = np.zeros(nt)
    z_end = np.zeros(5)
    loss_time = _api.simple_trace_orbit(state, nt, trajectory, times, z_end)
    result: dict[str, np.ndarray | float] = {
        "final_position": z_end,
        "loss_time": float(loss_time),
    }
    if return_trajectory:
        result["trajectory"] = np.ascontiguousarray(trajectory)
        result["times"] = times
    return result


def compute_canonical_frequencies(
    position: np.ndarray,
    *,
    integrator: str | int | None = None,
    n_periods: int = 1,
    max_steps: int = 10_000_000,
) -> dict[str, float | int | str]:
    """Bounce/transit period and mean toroidal precession of one orbit.

    A trapped cycle runs between same-direction ``v_par = 0`` crossings, a
    passing cycle between signed ``2*pi`` poloidal advances. With
    ``n_periods > 1`` the result is the mean and spread over several cycles.
    Periods are in seconds, angles in radians, frequencies in rad/s. Check
    ``status`` against the ``FREQ_*`` constants before using the values.
    """
    _require_init()
    state = _state(position)
    if n_periods < 1 or max_steps < 1:
        raise ValueError("n_periods and max_steps must be at least one")
    _use_integrator(integrator)
    values = np.zeros(6)
    metadata = np.zeros(5, dtype=np.int32)
    _api.simple_orbit_frequencies(state, int(n_periods), int(max_steps), values,
                                  metadata)
    return {
        "period": float(values[0]),
        "period_std": float(values[1]),
        "delta_phi": float(values[2]),
        "delta_phi_std": float(values[3]),
        "omega_b": float(values[4]),
        "omega_phi": float(values[5]),
        "status": int(metadata[0]),
        "orbit_class": _ORBIT_CLASS_NAMES.get(int(metadata[1]), "unknown"),
        "orbit_class_code": int(metadata[1]),
        "parallel_direction": int(metadata[2]),
        "n_periods": int(metadata[3]),
        "n_steps": int(metadata[4]),
    }


def invariants_from_state(position: np.ndarray) -> dict[str, float | int]:
    """Invariants ``(H0, J_perp, p_phi)`` of one state.

    ``H0 = (v/v0)**2`` and ``J_perp = H0 * (1 - lambda**2) / B``. ``p_phi`` is
    ``(c/q) P_phi`` in magnetic-flux units, zero on axis and positive outward
    in its flux term. It is not the radial label ``s``.
    """
    _require_init()
    values = np.zeros(6)
    status = _api.simple_invariants_from_state(_state(position), values)
    return {
        "h0": float(values[0]),
        "j_perp": float(values[1]),
        "p_phi": float(values[2]),
        "axis_flux_native": float(values[3]),
        "outward_sign": float(values[4]),
        "edge_flux_native": float(values[5]),
        "status": int(status),
    }


def states_from_invariants(
    h0: float, j_perp: float, p_phi: float, *, max_solutions: int = 64
) -> dict[str, np.ndarray | int]:
    """All start states on the axisymmetric cut with the given invariants.

    The cut is ``grad(psi) x grad(B) = 0``. States are ordered by major
    radius from the high- to the low-field side, as POTATO's ``R_c``.
    """
    _require_init()
    if not np.isfinite([h0, j_perp, p_phi]).all():
        raise ValueError("invariants must be finite")
    if h0 <= 0.0 or j_perp < 0.0 or max_solutions < 1:
        raise ValueError("need h0 > 0, j_perp >= 0 and max_solutions >= 1")
    states = np.zeros((5, max_solutions), order="F")
    metadata = np.zeros((3, max_solutions), dtype=np.int32, order="F")
    residuals = np.zeros(max_solutions)
    cylindrical = np.zeros((3, max_solutions), order="F")
    n_solutions, status = _api.simple_states_from_invariants(
        np.array([h0, j_perp, p_phi], dtype=np.float64), int(max_solutions),
        states, metadata, residuals, cylindrical,
    )
    n = int(n_solutions)
    kinds = metadata[2, :n].copy()
    return {
        "states": np.ascontiguousarray(states[:, :n].T),
        "cylindrical": np.ascontiguousarray(cylindrical[:, :n].T),
        "sigma": metadata[0, :n].copy(),
        "section_branch": metadata[1, :n].copy(),
        "section_kind_code": kinds,
        "section_kind": np.array([_SECTION_NAMES.get(int(k), "unknown")
                                  for k in kinds]),
        "residual": residuals[:n].copy(),
        "n_solutions": n,
        "status": int(status),
    }


def states_from_potato_invariants(
    h0: float,
    j_perp: float,
    psi_star: float,
    *,
    psi_axis: float,
    psi_edge: float,
    v0_ratio: float = 1.0,
    max_solutions: int = 64,
) -> dict[str, np.ndarray | int]:
    """Convert POTATO invariants and return the matching SIMPLE start states.

    ``v0_ratio`` is ``v0_POTATO / v0_SIMPLE``. ``psi_axis`` and ``psi_edge``
    use the units and gauge of POTATO's ``psi_star``.
    """
    if not (np.isfinite(v0_ratio) and v0_ratio > 0.0):
        raise ValueError("v0_ratio must be finite and positive")
    if not np.isfinite([psi_axis, psi_edge]).all() or psi_edge == psi_axis:
        raise ValueError("psi_axis and psi_edge must be finite and differ")
    converted = np.zeros(3)
    _api.simple_potato_invariants(
        np.array([h0, j_perp, psi_star], dtype=np.float64), float(psi_axis),
        float(psi_edge), float(v0_ratio), converted,
    )
    result = states_from_invariants(*converted, max_solutions=max_solutions)
    result["simple_invariants"] = converted
    return result


def trace_to_cut(
    position: np.ndarray,
    *,
    cut: str | int = "tip",
    integrator: str | int | None = None,
    max_events: int = 100,
) -> dict[str, np.ndarray | float | int | str]:
    """Trace one orbit to its next ``"tip"`` (``v_par = 0``), ``"toroidal"``
    field-period cut or ``"any"`` cut. Needs a symplectic integrator."""
    _require_init()
    cut_code = _CUTS.get(cut.lower()) if isinstance(cut, str) else int(cut)
    if cut_code not in _CUTS.values():
        raise ValueError(f"Unknown cut: {cut}")
    _use_integrator(integrator)
    if int(params.integmode) <= RK45:
        raise ValueError("trace_to_cut requires a symplectic integrator")
    if max_events < 1:
        raise ValueError("max_events must be at least one")
    cut_state = np.zeros(6)
    cut_type, status = _api.simple_trace_to_cut(
        _state(position), cut_code, int(max_events), cut_state
    )
    return {
        "state": cut_state[:5].copy(),
        "parallel_invariant": float(cut_state[5]),
        "cut_type": _CUT_NAMES.get(int(cut_type), "unknown"),
        "cut_type_code": int(cut_type),
        "status": int(status),
    }


__all__ = [
    "init", "sample", "sample_surface", "sample_volume", "load_particles",
    "trace_parallel", "classify_parallel", "classify_fast", "classify_fractal",
    "write_output", "trace_orbit", "compute_canonical_frequencies",
    "invariants_from_state", "states_from_invariants",
    "states_from_potato_invariants", "trace_to_cut", "params",
    "RK45", "EXPL_IMPL_EULER", "IMPL_EXPL_EULER", "MIDPOINT", "GAUSS1",
    "GAUSS2", "GAUSS3", "GAUSS4", "LOBATTO3",
    "ORBIT_CLASS_UNKNOWN", "ORBIT_CLASS_TRAPPED", "ORBIT_CLASS_PASSING",
    "FREQ_SUCCESS", "FREQ_INVALID_INPUT", "FREQ_ORBIT_LOST",
    "FREQ_INTEGRATOR_ERROR", "FREQ_MAX_STEPS",
    "INVARIANT_SUCCESS", "INVARIANT_INVALID_INPUT", "INVARIANT_NOT_AXISYMMETRIC",
    "INVARIANT_NO_INTERSECTION", "INVARIANT_ROOT_FAILURE", "INVARIANT_CAPACITY",
    "SECTION_EXTREMUM_UNKNOWN", "SECTION_B_MINIMUM", "SECTION_B_MAXIMUM",
    "CUT_ANY", "CUT_TIP", "CUT_TOROIDAL",
]

__all__ += ["cylindrical_to_reference", "reference_to_cylindrical"]
