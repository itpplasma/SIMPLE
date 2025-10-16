"""
Simulation result containers for the SIMPLE Python API.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from ._backend import SimulationArrays


@dataclass(slots=True)
class BatchResults:
    """
    Immutable view over SIMPLE simulation outputs.
    """

    _final_positions: np.ndarray
    _loss_times: np.ndarray
    _trap_parameter: Optional[np.ndarray]
    _perpendicular_invariant: Optional[np.ndarray]
    _tmax: float

    @classmethod
    def from_backend(cls, arrays: SimulationArrays) -> "BatchResults":
        return cls(
            _final_positions=arrays.final_positions,
            _loss_times=arrays.loss_times,
            _trap_parameter=arrays.trap_parameter,
            _perpendicular_invariant=arrays.perpendicular_invariant,
            _tmax=arrays.tmax,
        )

    @property
    def n_particles(self) -> int:
        return self._final_positions.shape[1]

    @property
    def tmax(self) -> float:
        return self._tmax

    @property
    def final_positions(self) -> np.ndarray:
        return self._final_positions

    @property
    def loss_times(self) -> np.ndarray:
        return self._loss_times

    @property
    def trap_parameter(self) -> Optional[np.ndarray]:
        return self._trap_parameter

    @property
    def perpendicular_invariant(self) -> Optional[np.ndarray]:
        return self._perpendicular_invariant

    def confined_mask(self, t_threshold: Optional[float] = None) -> np.ndarray:
        threshold = self._tmax if t_threshold is None else float(t_threshold)
        return self._loss_times >= threshold

    def lost_mask(self, t_threshold: Optional[float] = None) -> np.ndarray:
        return ~self.confined_mask(t_threshold)

    def confined(self, t_threshold: Optional[float] = None) -> np.ndarray:
        mask = self.confined_mask(t_threshold)
        return self._final_positions[:, mask]

    def lost(self, t_threshold: Optional[float] = None) -> dict[str, np.ndarray]:
        mask = self.lost_mask(t_threshold)
        return {
            "positions": self._final_positions[:, mask],
            "loss_times": self._loss_times[mask],
        }

    def to_dict(self) -> dict[str, np.ndarray | float]:
        return {
            "final_positions": np.array(self._final_positions, copy=True),
            "loss_times": np.array(self._loss_times, copy=True),
            "trap_parameter": None
            if self._trap_parameter is None
            else np.array(self._trap_parameter, copy=True),
            "perpendicular_invariant": None
            if self._perpendicular_invariant is None
            else np.array(self._perpendicular_invariant, copy=True),
            "tmax": self._tmax,
        }

    def summary(self) -> str:  # pragma: no cover - simple diagnostic formatting
        confined = self.confined_mask().sum()
        lost = self.lost_mask().sum()
        return (
            f"BatchResults(n_particles={self.n_particles}, "
            f"confined={confined}, lost={lost}, tmax={self._tmax})"
        )

    def __repr__(self) -> str:  # pragma: no cover - simple diagnostic formatting
        return self.summary()

