"""
Public entry-point for the cleaned SIMPLE Python API.
"""

from __future__ import annotations

import os
import urllib.request
from contextlib import nullcontext
from dataclasses import dataclass
from enum import IntEnum
from pathlib import Path
from typing import Any, Dict, Optional

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

DEFAULT_VMEC_URL = (
    "https://github.com/hiddenSymmetries/simsopt/raw/master/"
    "tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
)


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


def ensure_example_vmec(
    destination: str | Path = "wout.nc", *, source_url: str = DEFAULT_VMEC_URL
) -> Path:
    """Download the reference VMEC file if it is not yet available."""

    dest = Path(destination).expanduser()
    if dest.is_dir():
        dest = dest / "wout.nc"

    if dest.exists():
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(source_url, dest)
    return dest


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

    _backend.run_simulation(batch.positions, tmax, integrator_code, verbose=verbose)
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


class JParallelClass(IntEnum):
    """Classification labels for the :math:`J_\\parallel` heuristic."""

    UNCLASSIFIED = 0
    REGULAR = 1
    STOCHASTIC = 2


class TopologyClass(IntEnum):
    """Classification labels for the topological ideal/non-ideal heuristic."""

    UNCLASSIFIED = 0
    IDEAL = 1
    NON_IDEAL = 2


class MinkowskiClass(IntEnum):
    """Classification labels for the Minkowski fractal-dimension heuristic."""

    UNCLASSIFIED = 0
    REGULAR = 1
    STOCHASTIC = 2


@dataclass(slots=True)
class ClassificationResult:
    """Container for classifier outputs."""

    j_parallel: np.ndarray
    topology: np.ndarray
    minkowski: np.ndarray
    trap_parameter: np.ndarray
    loss_times: np.ndarray

    def as_enums(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the classification arrays cast to their enum types."""
        return (
            self.j_parallel.astype(JParallelClass),
            self.topology.astype(TopologyClass),
            self.minkowski.astype(MinkowskiClass),
        )

    def counts(self) -> Dict[str, Dict[int, int]]:
        """Return histogram-style counts for each classifier."""

        def hist(values: np.ndarray) -> Dict[int, int]:
            labels, counts = np.unique(values.astype(int), return_counts=True)
            return {int(label): int(count) for label, count in zip(labels, counts)}

        return {
            "j_parallel": hist(self.j_parallel),
            "topology": hist(self.topology),
            "minkowski": hist(self.minkowski),
        }


# Backwards compatibility name
FastClassificationResult = ClassificationResult


def classify_fast(
    particles: ParticleBatch | np.ndarray,
    *,
    tcut: float = 0.1,
    vmec_file: str | Path | None = None,
    integrator: str | int | None = None,
    write_files: bool = False,
    include_minkowski: bool = True,
) -> ClassificationResult:
    """Run the fast classifiers using an ephemeral :class:`SimpleSession`."""

    vmec = vmec_file or current_vmec()
    if vmec is None:
        raise RuntimeError("VMEC equilibrium not initialized. Provide vmec_file.")

    session = SimpleSession(vmec, integrator=integrator)
    return session.classify_fast(
        particles,
        classification_time=tcut,
        legacy_files=write_files,
        output_dir=None,
        assume_passing_confined=True,
        integrator=integrator,
        include_minkowski=include_minkowski,
    )


class SimpleSession:
    """High-level helper wrapping VMEC loading and API calls."""

    def __init__(
        self,
        vmec_file: str | Path,
        *,
        integrator: str | int | None = None,
    ) -> None:
        self.vmec_file = Path(vmec_file).expanduser()
        load_vmec(self.vmec_file)
        self._integrator = integrator

    def _integrator_choice(self, integrator: str | int | None) -> str | int | None:
        return integrator if integrator is not None else self._integrator

    def sample_surface(
        self, n_particles: int, *, surface: float = DEFAULT_SURFACE
    ) -> ParticleBatch:
        batch = ParticleBatch(n_particles)
        batch.initialize_from_samplers(self.vmec_file, method="surface", s=surface)
        return batch

    def sample_volume(
        self,
        n_particles: int,
        *,
        s_inner: float = DEFAULT_S_INNER,
        s_outer: float = DEFAULT_S_OUTER,
    ) -> ParticleBatch:
        batch = ParticleBatch(n_particles)
        batch.initialize_from_samplers(
            self.vmec_file, method="volume", s_inner=s_inner, s_outer=s_outer
        )
        return batch

    def load_particles(self, particle_file: str | Path) -> ParticleBatch:
        data = load_particle_file(self.vmec_file, particle_file)
        return ParticleBatch.from_fortran_arrays(data)

    def trace(
        self,
        particles: ParticleBatch | np.ndarray,
        *,
        tmax: float = DEFAULT_TMAX,
        integrator: str | int | None = None,
        verbose: bool = False,
    ) -> BatchResults:
        return trace_orbits(
            particles,
            tmax=tmax,
            integrator=self._integrator_choice(integrator),
            vmec_file=self.vmec_file,
            verbose=verbose,
        )

    def classify_fast(
        self,
        particles: ParticleBatch | np.ndarray,
        *,
        classification_time: float = 0.1,
        integrator: str | int | None = None,
        assume_passing_confined: bool = True,
        legacy_files: bool = False,
        output_dir: Path | None = None,
        include_minkowski: bool = True,
    ) -> ClassificationResult:
        """
        Run the fast classifiers and return their labels.

        Parameters
        ----------
        particles:
            Particle batch or array containing initial conditions.
        classification_time:
            Duration of the fast classification trace.
        integrator:
            Optional integrator override.
        assume_passing_confined:
            If ``True`` (default) passing particles are marked as confined without
            tracing the full orbit.
        legacy_files:
            When ``True`` the legacy ``fort.*`` outputs are produced in ``output_dir``.
        output_dir:
            Target directory for legacy outputs.
        include_minkowski:
            Set to ``False`` to skip the Minkowski fractal-dimension classifier and run
            only the fast J_parallel/topology heuristics.
        """
        batch = _as_batch(particles)
        n_particles = batch.n_particles

        overrides: Dict[str, float | int | bool] = {
            "fast_class": True,
            "class_plot": legacy_files,
        }
        if include_minkowski:
            overrides["tcut"] = classification_time
        else:
            overrides["tcut"] = -1.0

        if assume_passing_confined:
            overrides["notrace_passing"] = 1

        trace_params = get_parameters("trace_time")
        if trace_params["trace_time"] < classification_time:
            overrides["trace_time"] = classification_time

        with _backend.working_directory(output_dir):
            with _backend.temporary_parameters(**overrides):
                trace_orbits(
                    batch,
                    tmax=classification_time,
                    integrator=self._integrator_choice(integrator),
                    vmec_file=None,
                    verbose=False,
                )
                iclass = _backend.snapshot_classification(n_particles)
                trap_parameter = _backend.snapshot_trap_parameter(n_particles)
                loss_times = _backend.snapshot_loss_times(n_particles)

        if not include_minkowski:
            iclass[2, :] = 0

        return ClassificationResult(
            j_parallel=iclass[0, :].astype(np.int64, copy=False),
            topology=iclass[1, :].astype(np.int64, copy=False),
            minkowski=iclass[2, :].astype(np.int64, copy=False),
            trap_parameter=trap_parameter,
            loss_times=loss_times,
        )


temporary_parameters = _backend.temporary_parameters


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
    "JParallelClass",
    "TopologyClass",
    "MinkowskiClass",
    "ClassificationResult",
    "FastClassificationResult",
    "SimpleSession",
    "ensure_example_vmec",
    "temporary_parameters",
    "classify_fast",
    "load_particle_file",
]
