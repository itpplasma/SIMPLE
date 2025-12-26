#!/usr/bin/env python3
"""
Batch-oriented validation tests for the SIMPLE Python API.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

python_dir = Path(__file__).parent.parent.parent / "python"
sys.path.insert(0, str(python_dir))

sys.modules.pop("pysimple", None)
pysimple = pytest.importorskip("pysimple", reason="pysimple module not available")


class TestSoALayout:
    def test_numpy_layout(self, vmec_file: str):
        # Explicitly set ntestpart to avoid interference
        pysimple.init(vmec_file, deterministic=True, ntestpart=32)
        positions = pysimple.sample_surface(32, s=0.55)

        assert positions.flags.c_contiguous
        assert positions.shape == (5, 32)

        s_coords = positions[0, :]
        phi_coords = positions[2, :]
        assert s_coords.size == 32
        assert phi_coords.size == 32

    def test_zero_copy_write(self, vmec_file: str):
        # Explicitly set ntestpart to avoid interference
        pysimple.init(vmec_file, deterministic=True, ntestpart=4)
        positions = pysimple.sample_surface(4, s=0.2)

        original_value = positions[1, 0]
        positions[1, 0] += 0.25
        assert positions[1, 0] == original_value + 0.25


class TestBatchResults:
    def test_confined_statistics(self, vmec_file: str):
        # Explicitly set ntestpart to avoid interference
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=10)
        particles = pysimple.sample_surface(10, s=0.5)
        results = pysimple.trace_parallel(particles, integrator="midpoint")

        confined_mask = results['loss_times'] >= 1e-4
        assert confined_mask.shape == (10,)
        assert confined_mask.dtype == bool

        lost_mask = results['loss_times'] < 1e-4
        lost_positions = results['final_positions'][:, lost_mask]
        assert lost_positions.shape[0] == 5

    def test_to_dict(self, vmec_file: str):
        # Explicitly set ntestpart to avoid interference
        pysimple.init(vmec_file, deterministic=True, trace_time=5e-5, ntestpart=3)
        particles = pysimple.sample_surface(3, s=0.3)
        results = pysimple.trace_parallel(particles, integrator="midpoint")

        assert isinstance(results, dict)
        assert set(results.keys()) == {
            "final_positions",
            "loss_times",
            "trap_parameter",
            "perpendicular_invariant",
        }
        assert results["final_positions"].shape == (5, 3)


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
