"""
Public entry-point for the cleaned SIMPLE Python API.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import numpy as np

from . import _backend
from .particles import ParticleBatch
from .results import BatchResults
from .samplers import SurfaceSampler, VolumeSampler, load_particle_file

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

DEFAULT_TMAX = 0.1
DEFAULT_SURFACE = 0.5
DEFAULT_S_INNER = 0.1
DEFAULT_S_OUTER = 0.9
DEFAULT_INTEGRATOR = MIDPOINT

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


def load_vmec(
    vmec_file: str | Path,
    *,
    ns_s: int = _backend.DEFAULT_NS_S,
    ns_tp: int = _backend.DEFAULT_NS_TP,
    multharm: int = _backend.DEFAULT_MULTHARM,
    integrator: int | None = None,
    force: bool = False,
) -> None:
    """Load or re-load a VMEC equilibrium into the Fortran backend."""
    _backend.ensure_vmec_loaded(
        vmec_file,
        ns_s=ns_s,
        ns_tp=ns_tp,
        multharm=multharm,
        integrator=integrator,
        force=force,
    )


def current_vmec() -> Optional[str]:
    """Return the currently loaded VMEC equilibrium path, if any."""
    return _backend.current_vmec()


def _resolve_integrator(integrator: str | int | None) -> int:
    if integrator is None:
        return DEFAULT_INTEGRATOR
    if isinstance(integrator, int):
        return integrator
    key = integrator.lower()
    if key not in _INTEGRATOR_ALIASES:
        raise ValueError(f"Unknown integrator '{integrator}'")
    return _INTEGRATOR_ALIASES[key]


def _as_batch(particles: ParticleBatch | np.ndarray) -> ParticleBatch:
    if isinstance(particles, ParticleBatch):
        return particles
    return ParticleBatch.from_array(np.asarray(particles))


def trace_orbits(
    particles: ParticleBatch | np.ndarray,
    *,
    tmax: float = DEFAULT_TMAX,
    integrator: str | int | None = None,
    vmec_file: str | Path | None = None,
    verbose: bool = False,
) -> BatchResults:
    """
    Trace particle orbits using the SIMPLE Fortran backend.
    """
    if vmec_file is not None:
        load_vmec(vmec_file)
    else:
        _backend.assert_vmec_loaded()

    batch = _as_batch(particles)
    integrator_code = _resolve_integrator(integrator)

    _backend.update_start_positions(batch.positions)
    _backend.run_simulation(tmax, integrator_code, verbose=verbose)
    arrays = _backend.collect_results(tmax)

    return BatchResults.from_backend(arrays)


def set_parameters(**kwargs: Any) -> None:
    """Set Fortran parameters by delegating to the backend module."""
    _backend.set_params(**kwargs)


def get_parameters(*names: str) -> dict[str, Any]:
    """Fetch one or more Fortran parameters."""
    return _backend.get_params(*names)


def get_confined(results: BatchResults, t_threshold: float | None = None) -> np.ndarray:
    """Convenience wrapper for :meth:`BatchResults.confined`."""
    return results.confined(t_threshold)


def get_lost(
    results: BatchResults, t_threshold: float | None = None
) -> dict[str, np.ndarray]:
    """Convenience wrapper for :meth:`BatchResults.lost`."""
    return results.lost(t_threshold)


__all__ = [
    "ParticleBatch",
    "BatchResults",
    "SurfaceSampler",
    "VolumeSampler",
    "trace_orbits",
    "load_vmec",
    "current_vmec",
    "set_parameters",
    "get_parameters",
    "get_confined",
    "get_lost",
    "DEFAULT_TMAX",
    "DEFAULT_SURFACE",
    "DEFAULT_S_INNER",
    "DEFAULT_S_OUTER",
    "DEFAULT_INTEGRATOR",
    "RK45",
    "EXPL_IMPL_EULER",
    "IMPL_EXPL_EULER",
    "MIDPOINT",
    "GAUSS1",
    "GAUSS2",
    "GAUSS3",
    "GAUSS4",
    "LOBATTO3",
    "load_particle_file",
]
