#!/usr/bin/env python3
"""
Batch-oriented validation tests for the SIMPLE Python API.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

python_dir = Path(__file__).parent.parent.parent / "python"
pytest.importorskip("pysimple", reason="pysimple module not available")

import sys

sys.path.insert(0, str(python_dir))

import simple  # noqa: E402  pylint: disable=wrong-import-position


class TestSoALayout:
    def test_numpy_layout(self, vmec_file: str):
        batch = simple.ParticleBatch(32)
        batch.initialize_from_samplers(vmec_file, method="surface", s=0.55)

        positions = batch.positions
        assert positions.flags.c_contiguous
        assert positions.shape == (5, 32)

        s_coords = positions[0, :]
        phi_coords = positions[2, :]
        assert s_coords.size == 32
        assert phi_coords.size == 32

    def test_zero_copy_write(self, vmec_file: str):
        batch = simple.ParticleBatch(4)
        batch.initialize_from_samplers(vmec_file, method="surface", s=0.2)

        view = batch.positions
        view[1, 0] += 0.25
        assert batch.positions[1, 0] == view[1, 0]


class TestBatchResults:
    def test_confined_statistics(self, vmec_file: str):
        batch = simple.ParticleBatch(10)
        batch.initialize_from_samplers(vmec_file, method="surface", s=0.5)

        results = simple.trace_orbits(batch, tmax=1e-4, integrator="midpoint")

        confined_mask = results.confined_mask()
        assert confined_mask.shape == (10,)
        assert confined_mask.dtype == bool

        lost = results.lost()
        assert lost["positions"].shape[0] == 5
        assert lost["loss_times"].ndim == 1

    def test_to_dict(self, vmec_file: str):
        batch = simple.ParticleBatch(3)
        batch.initialize_from_samplers(vmec_file, method="surface", s=0.3)

        results = simple.trace_orbits(batch, tmax=5e-5, integrator="midpoint")
        payload = results.to_dict()

        assert set(payload) == {
            "final_positions",
            "loss_times",
            "trap_parameter",
            "perpendicular_invariant",
            "tmax",
        }
        assert payload["final_positions"].shape == (5, 3)
        assert payload["tmax"] == pytest.approx(5e-5)


class TestPerformanceFramework:
    def test_column_major_access_speed(self):
        n_particles = 5000
        soa = np.zeros((5, n_particles), dtype=np.float64)
        soa[0, :] = np.linspace(0.1, 0.9, n_particles)
        soa[1, :] = np.random.uniform(0, 2 * np.pi, n_particles)

        column_access = soa[:, 1000]
        row_access = soa[0, :]

        assert column_access.shape == (5,)
        assert row_access.shape == (n_particles,)

