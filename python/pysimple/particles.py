"""
Particle batch abstractions for the cleaned SIMPLE Python API.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from . import _backend


def _coerce_positions(positions: np.ndarray) -> np.ndarray:
    """
    Ensure positions follow the canonical (5, n_particles) SoA layout.
    """
    array = np.asarray(positions, dtype=np.float64)
    if array.ndim != 2:
        raise ValueError("Particle arrays must be two-dimensional (5, n_particles)")

    if array.shape[0] == 5:
        return np.ascontiguousarray(array, dtype=np.float64)

    if array.shape[1] == 5:
        return np.ascontiguousarray(array.T, dtype=np.float64)

    raise ValueError(
        f"Particle arrays must have shape (5, N) or (N, 5); received {array.shape}"
    )


@dataclass(slots=True)
class ParticleBatch:
    """
    Structure-of-arrays container for particle phase-space coordinates.

    The batch stores data in the canonical SIMPLE order:
    ``(s, theta, phi, v_parallel, magnetic_moment)``.
    """

    _positions: np.ndarray

    def __init__(self, n_particles: int, *, dtype=np.float64) -> None:
        n_particles = int(n_particles)
        if n_particles < 0:
            raise ValueError("n_particles must be non-negative")

        self._positions = np.zeros((5, n_particles), dtype=dtype, order="C")

    @property
    def positions(self) -> np.ndarray:
        """Direct, writable view of the particle coordinates."""
        return self._positions

    @property
    def n_particles(self) -> int:
        """The number of particles in the batch."""
        return self._positions.shape[1]

    def copy(self) -> "ParticleBatch":
        """Deep copy of the particle batch."""
        clone = ParticleBatch(self.n_particles)
        clone._positions[...] = self._positions
        return clone

    def to_numpy(self, *, copy: bool = True) -> np.ndarray:
        """Return the positions as a NumPy array."""
        if copy:
            return np.array(self._positions, copy=True)
        return self._positions

    def initialize_from_samplers(
        self,
        vmec_file: Optional[str | Path] = None,
        *,
        method: str = "surface",
        s: float = 0.5,
        s_inner: float = 0.1,
        s_outer: float = 0.9,
        particle_file: Optional[str | Path] = None,
    ) -> "ParticleBatch":
        """
        Populate the batch using the Fortran sampling routines.
        """
        if vmec_file is not None:
            _backend.ensure_vmec_loaded(vmec_file)
        else:
            _backend.assert_vmec_loaded()

        if method == "surface":
            sampled = _backend.surface_sample(self.n_particles, s)
        elif method == "volume":
            sampled = _backend.volume_sample(self.n_particles, s_inner, s_outer)
        elif method == "file":
            if particle_file is None:
                raise ValueError("particle_file must be provided when method='file'")
            sampled = _backend.load_particle_file(particle_file)
            if sampled.shape[1] != self.n_particles:
                self._positions = np.zeros(
                    (5, sampled.shape[1]), dtype=np.float64, order="C"
                )
        else:
            raise ValueError(f"Unsupported sampling method '{method}'")

        if sampled.shape != self._positions.shape:
            self._positions = np.ascontiguousarray(sampled, dtype=np.float64)
        else:
            self._positions[...] = sampled

        return self

    @classmethod
    def from_array(cls, positions: np.ndarray) -> "ParticleBatch":
        """Construct a batch from a raw numpy array."""
        coerced = _coerce_positions(positions)
        batch = cls(coerced.shape[1])
        batch._positions[...] = coerced
        return batch

    @classmethod
    def from_fortran_arrays(cls, positions: np.ndarray) -> "ParticleBatch":
        """
        Construct a batch from arrays returned by the Fortran samplers.
        """
        return cls.from_array(positions)

    # NumPy interoperability -------------------------------------------------
    def __array__(self, dtype=None) -> np.ndarray:
        return np.asarray(self._positions, dtype=dtype)

    def __repr__(self) -> str:  # pragma: no cover - trivial
        return f"ParticleBatch(n_particles={self.n_particles})"
