#!/usr/bin/env python3
"""
Tests for the cleaned SIMPLE Python API.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

python_dir = Path(__file__).parent.parent.parent / "python"
pytest.importorskip("pysimple", reason="pysimple module not available")

import sys

sys.path.insert(0, str(python_dir))

examples_dir = Path(__file__).resolve().parents[2] / "examples"
sys.path.insert(0, str(examples_dir))

import pysimple  # noqa: E402  pylint: disable=wrong-import-position


class TestIntegratorSelection:
    def test_alias_resolution(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True, trace_time=5e-5)
        particles = pysimple.sample_surface(2, s=0.42)

        pysimple.trace_parallel(particles.copy(), integrator="symplectic_midpoint")
        pysimple.trace_parallel(particles.copy(), integrator="rk45")

        with pytest.raises(ValueError):
            pysimple.trace_parallel(particles.copy(), integrator="does_not_exist")


class TestSurfaceSampler:
    def test_surface_sampler_shape(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True)
        particles = pysimple.sample_surface(8, s=0.4)
        assert particles.shape == (5, 8)
        assert particles.flags.c_contiguous


class TestSampleSurface:
    def test_sample_surface(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True)
        particles = pysimple.sample_surface(12, s=0.35)

        assert particles.shape == (5, 12)
        assert not np.allclose(particles, 0.0)

        original = particles[0, 0]
        particles[0, 0] = original + 1.0
        assert particles[0, 0] == original + 1.0

    def test_particles_can_be_copied(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True)
        particles = pysimple.sample_surface(6, s=0.5)

        clone = particles.copy()
        np.testing.assert_allclose(particles, clone)
        clone[0, 0] += 0.1
        assert clone[0, 0] != particles[0, 0]


class TestTraceOrbits:
    def test_trace_returns_dict_results(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4)
        particles = pysimple.sample_surface(4, s=0.3)
        results = pysimple.trace_parallel(particles)

        assert isinstance(results, dict)
        assert results['final_positions'].shape == (5, 4)
        assert results['loss_times'].shape == (4,)

    def test_confined_and_lost_filters(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True, trace_time=5e-5)
        particles = pysimple.sample_surface(3, s=0.3)
        results = pysimple.trace_parallel(particles)

        confined_mask = results['loss_times'] >= 5e-5
        lost_mask = results['loss_times'] < 5e-5

        assert confined_mask.shape == (3,)
        assert lost_mask.shape == (3,)
