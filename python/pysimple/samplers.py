"""
High level helpers around the SIMPLE Fortran samplers.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from . import _backend


@dataclass(slots=True)
class SurfaceSampler:
    """
    Wrapper for surface-based particle initialisation.
    """

    vmec_file: str
    ns_s: int = _backend.DEFAULT_NS_S
    ns_tp: int = _backend.DEFAULT_NS_TP
    multharm: int = _backend.DEFAULT_MULTHARM

    def __post_init__(self) -> None:
        _backend.ensure_vmec_loaded(
            self.vmec_file, ns_s=self.ns_s, ns_tp=self.ns_tp, multharm=self.multharm
        )

    def sample_surface_fieldline(
        self, n_particles: int, *, s: float = 0.5
    ) -> np.ndarray:
        return _backend.surface_sample(int(n_particles), float(s))

    def initialize_batch(
        self, batch: "ParticleBatch", *, s: float = 0.5
    ) -> "ParticleBatch":
        from .particles import ParticleBatch  # Local import to avoid cycle

        if not isinstance(batch, ParticleBatch):
            raise TypeError("batch must be a ParticleBatch instance")
        batch.initialize_from_samplers(method="surface", s=s)
        return batch


@dataclass(slots=True)
class VolumeSampler:
    """
    Wrapper for volume-based particle sampling.
    """

    vmec_file: str
    ns_s: int = _backend.DEFAULT_NS_S
    ns_tp: int = _backend.DEFAULT_NS_TP
    multharm: int = _backend.DEFAULT_MULTHARM

    def __post_init__(self) -> None:
        _backend.ensure_vmec_loaded(
            self.vmec_file, ns_s=self.ns_s, ns_tp=self.ns_tp, multharm=self.multharm
        )

    def sample_volume(
        self, n_particles: int, *, s_inner: float = 0.1, s_outer: float = 0.9
    ) -> np.ndarray:
        return _backend.volume_sample(
            int(n_particles), float(s_inner), float(s_outer)
        )

    def initialize_batch(
        self,
        batch: "ParticleBatch",
        *,
        s_inner: float = 0.1,
        s_outer: float = 0.9,
    ) -> "ParticleBatch":
        from .particles import ParticleBatch

        if not isinstance(batch, ParticleBatch):
            raise TypeError("batch must be a ParticleBatch instance")
        batch.initialize_from_samplers(
            method="volume", s_inner=s_inner, s_outer=s_outer
        )
        return batch


def load_particle_file(
    vmec_file: Optional[str | Path],
    particle_file: str | Path,
) -> np.ndarray:
    """
    Load particles from disk after ensuring VMEC is ready.
    """
    if vmec_file is not None:
        _backend.ensure_vmec_loaded(vmec_file)
    else:
        _backend.assert_vmec_loaded()
    return _backend.load_particle_file(particle_file)

