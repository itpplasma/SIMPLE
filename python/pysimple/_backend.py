"""
Internal helper utilities for the cleaned SIMPLE Python API.

This module centralizes all direct interactions with the ``simple_backend`` f90wrap
bindings so the public API can remain small and well tested.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

import os

import numpy as np

try:
    import simple_backend as pysimple  # type: ignore
except ImportError as exc:  # pragma: no cover - exercised in test skips
    raise ImportError(
        "simple_backend module not found. Build SIMPLE with Python bindings enabled."
    ) from exc


# Default spline orders used by the Fortran driver (examples/simple_full.in)
DEFAULT_NS_S = 5
DEFAULT_NS_TP = 5
DEFAULT_MULTHARM = 5


@dataclass(slots=True)
class SimulationArrays:
    """Container for numpy views copied out of the Fortran arrays."""

    final_positions: np.ndarray
    loss_times: np.ndarray
    trap_parameter: Optional[np.ndarray]
    perpendicular_invariant: Optional[np.ndarray]
    tmax: float

    @property
    def n_particles(self) -> int:
        return self.final_positions.shape[1]


_current_vmec: Optional[str] = None
_SAMPLERS = pysimple.Samplers()


def current_vmec() -> Optional[str]:
    """Return the currently loaded VMEC file, if any."""
    return _current_vmec


def ensure_vmec_loaded(
    vmec_file: str | Path,
    *,
    ns_s: int = DEFAULT_NS_S,
    ns_tp: int = DEFAULT_NS_TP,
    multharm: int = DEFAULT_MULTHARM,
    integrator: Optional[int] = None,
    force: bool = False,
) -> None:
    """
    Initialize SIMPLE for a VMEC equilibrium if not already loaded.

    Parameters
    ----------
    vmec_file:
        Path to the VMEC NetCDF equilibrium.
    ns_s, ns_tp, multharm:
        Spline order configuration passed to the Fortran initialization
        routine. Defaults mirror the canonical SIMPLE configuration.
    integrator:
        Optional integration mode to seed during initialization.
    force:
        Re-initialize even if the requested file matches the cached VMEC.
    """
    global _current_vmec

    vmec_path = str(Path(vmec_file).expanduser().resolve())

    if not force and _current_vmec == vmec_path:
        return

    if integrator is not None:
        pysimple.params.integmode = int(integrator)

    pysimple.params.netcdffile = vmec_path

    tracer = pysimple.simple.Tracer()
    pysimple.simple_main.init_field(
        tracer,
        vmec_path,
        int(ns_s),
        int(ns_tp),
        int(multharm),
        pysimple.params.integmode,
    )
    pysimple.params.params_init()

    _current_vmec = vmec_path


def get_field_type() -> int:
    """Return the currently active field-type flag."""
    return int(pysimple.velo_mod.isw_field_type)


def set_field_type(value: int, *, reload: bool = True) -> None:
    """
    Update the field-type flag, optionally forcing a VMEC reload.
    """
    target = int(value)
    if get_field_type() == target and not reload:
        return

    pysimple.velo_mod.isw_field_type = target

    if reload and _current_vmec is not None:
        ensure_vmec_loaded(_current_vmec, force=True)


def assert_vmec_loaded() -> None:
    """Raise if no VMEC equilibrium has been initialized."""
    if _current_vmec is None:
        raise RuntimeError("VMEC equilibrium not initialized. Call load_vmec first.")


def configure_batch(n_particles: int) -> None:
    """
    Resize Fortran-side arrays for the requested particle count.
    """
    assert_vmec_loaded()

    count = int(n_particles)
    pysimple.params.ntestpart = count
    pysimple.params.reallocate_arrays()


def collect_results(tmax: float) -> SimulationArrays:
    """
    Snapshot the Fortran output arrays after a simulation run.
    """
    assert_vmec_loaded()

    n_particles = int(pysimple.params.ntestpart)
    final_positions = np.ascontiguousarray(
        pysimple.params.zend[:, :n_particles], dtype=np.float64
    )
    loss_times = np.ascontiguousarray(
        pysimple.params.times_lost[:n_particles], dtype=np.float64
    )

    trap_parameter = None
    if hasattr(pysimple.params, "trap_par"):
        trap_parameter = np.ascontiguousarray(
            pysimple.params.trap_par[:n_particles], dtype=np.float64
        )

    perpendicular_invariant = None
    if hasattr(pysimple.params, "perp_inv"):
        perpendicular_invariant = np.ascontiguousarray(
            pysimple.params.perp_inv[:n_particles], dtype=np.float64
        )

    return SimulationArrays(
        final_positions=final_positions,
        loss_times=loss_times,
        trap_parameter=trap_parameter,
        perpendicular_invariant=perpendicular_invariant,
        tmax=float(tmax),
    )


def _callback_module_available() -> bool:
    return hasattr(pysimple, "callback") and hasattr(pysimple.callback, "set_output_orbits_macrostep")


def set_macrostep_output(enabled: bool) -> None:
    if not _callback_module_available():
        raise RuntimeError("Macrostep orbit output not available in this build of SIMPLE")
    pysimple.callback.set_output_orbits_macrostep(bool(enabled))
    if hasattr(pysimple.params, "output_orbits_macrostep"):
        pysimple.params.output_orbits_macrostep = bool(enabled)


def get_macrostep_output() -> bool:
    if not _callback_module_available():
        return False
    if hasattr(pysimple.params, "output_orbits_macrostep"):
        return bool(pysimple.params.output_orbits_macrostep)
    return bool(pysimple.callback.get_output_orbits_macrostep())


def run_simulation(
    positions: np.ndarray,
    tmax: float,
    integrator_code: int,
    verbose: bool = False,
) -> None:
    """
    Execute the SIMPLE Fortran integrator for the current batch.
    """
    assert_vmec_loaded()

    batch_positions = np.ascontiguousarray(positions, dtype=np.float64)
    n_particles = batch_positions.shape[1]

    configure_batch(n_particles)
    pysimple.params.trace_time = float(tmax)
    pysimple.params.integmode = int(integrator_code)
    pysimple.params.params_init()
    pysimple.params.zstart[:, :n_particles] = batch_positions

    tracer = pysimple.simple.Tracer()
    pysimple.simple_main.run(tracer)

    if verbose:
        pysimple.simple_main.print_parameters()


def surface_sample(n_particles: int, s: float) -> np.ndarray:
    """
    Generate particles on a flux surface using the Fortran sampler routines.
    """
    assert_vmec_loaded()
    configure_batch(n_particles)

    pysimple.params.startmode = 2
    pysimple.params.sbeg[0] = float(s)
    pysimple.params.generate_start_only = False

    _SAMPLERS.sample_surface_fieldline(pysimple.params.zstart)

    return np.ascontiguousarray(pysimple.params.zstart[:, :n_particles], dtype=np.float64)


def volume_sample(n_particles: int, s_inner: float, s_outer: float) -> np.ndarray:
    """
    Generate particles in a volume using the Fortran sampler routines.
    """
    assert_vmec_loaded()
    configure_batch(n_particles)

    pysimple.params.startmode = 3
    pysimple.params.generate_start_only = False

    _SAMPLERS.sample_volume_single(
        pysimple.params.zstart, float(s_inner), float(s_outer)
    )

    return np.ascontiguousarray(pysimple.params.zstart[:, :n_particles], dtype=np.float64)


def load_particle_file(particle_file: str | Path) -> np.ndarray:
    """
    Load particles from a text file using the Fortran sampler.
    """
    assert_vmec_loaded()

    particle_path = Path(particle_file).expanduser().resolve()

    with particle_path.open("r", encoding="utf-8") as handle:
        n_particles = sum(
            1 for line in handle if line.strip() and not line.lstrip().startswith("#")
        )

    if n_particles == 0:
        return np.zeros((5, 0), dtype=np.float64, order="C")

    configure_batch(n_particles)
    pysimple.params.startmode = 1
    pysimple.params.generate_start_only = False

    _SAMPLERS.sample_read(pysimple.params.zstart, str(particle_path))

    return np.ascontiguousarray(pysimple.params.zstart[:, :n_particles], dtype=np.float64)


def set_params(**kwargs: float | int | bool) -> None:
    """Thin wrapper around ``pysimple.params`` attribute assignment."""
    assert_vmec_loaded()

    for key, value in kwargs.items():
        if not hasattr(pysimple.params, key):
            raise ValueError(f"Unknown SIMPLE parameter '{key}'")
        setattr(pysimple.params, key, value)


def get_params(*names: str) -> dict[str, float | int | bool]:
    """Fetch parameters from the Fortran module."""
    assert_vmec_loaded()

    if not names:
        raise ValueError("Specify at least one parameter name")
    return {name: getattr(pysimple.params, name) for name in names}


def snapshot_classification(n_particles: int) -> np.ndarray:
    """Copy the ``iclass`` array for the leading particles."""
    assert_vmec_loaded()
    return np.array(
        pysimple.params.iclass[:, :n_particles],
        copy=True,
        dtype=np.int64,
    )


def snapshot_trap_parameter(n_particles: int) -> np.ndarray:
    """Copy the trapped-parameter array for the leading particles."""
    assert_vmec_loaded()
    return np.array(pysimple.params.trap_par[:n_particles], copy=True)


def snapshot_loss_times(n_particles: int) -> np.ndarray:
    """Copy the loss-time array for the leading particles."""
    assert_vmec_loaded()
    return np.array(pysimple.params.times_lost[:n_particles], copy=True)


def snapshot_start_positions(n_particles: int) -> np.ndarray:
    """Copy the current ``zstart`` array for the leading particles."""
    assert_vmec_loaded()
    return np.array(pysimple.params.zstart[:, :n_particles], copy=True)


@contextmanager
def temporary_parameters(**overrides: float | int | bool) -> Iterator[None]:
    """Temporarily modify Fortran parameters within a context."""
    if not overrides:
        yield
        return

    names = tuple(overrides.keys())
    original = get_params(*names)
    set_params(**overrides)
    try:
        yield
    finally:
        set_params(**original)


@contextmanager
def working_directory(path: Optional[Path]) -> Iterator[None]:
    """Context manager that changes the current working directory."""
    if path is None:
        yield
        return

    target = Path(path).expanduser().resolve()
    target.mkdir(parents=True, exist_ok=True)
    prev = Path.cwd()
    os.chdir(target)
    try:
        yield
    finally:
        os.chdir(prev)


@contextmanager
def macrostep_output(enabled: bool = True) -> Iterator[None]:
    if not _callback_module_available():
        if enabled:
            raise RuntimeError("Macrostep orbit output not available in this build of SIMPLE")
        yield
        return

    previous = get_macrostep_output()
    set_macrostep_output(enabled)
    try:
        yield
    finally:
        set_macrostep_output(previous)


@contextmanager
def field_type(value: int) -> Iterator[None]:
    """
    Temporarily switch the field-type flag in a safe manner.
    """
    previous = get_field_type()
    target = int(value)

    if previous == target:
        yield
        return

    set_field_type(target, reload=True)
    try:
        yield
    finally:
        set_field_type(previous, reload=True)
