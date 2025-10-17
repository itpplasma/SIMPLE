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

import simple  # noqa: E402  pylint: disable=wrong-import-position
from classify_fast import classify_fast_example  # noqa: E402
from simple_api import run_trace_example  # noqa: E402


@pytest.fixture(scope="module")
def sampler(vmec_file: str) -> simple.SurfaceSampler:
    return simple.SurfaceSampler(vmec_file)


class TestIntegratorSelection:
    def test_alias_resolution(self, vmec_file: str):
        batch = simple.ParticleBatch(2)
        batch.initialize_from_samplers(vmec_file, method="surface", s=0.42)

        simple.trace_orbits(batch.copy(), tmax=5e-5, integrator="symplectic_midpoint")
        simple.trace_orbits(batch.copy(), tmax=5e-5, integrator="rk45")

        with pytest.raises(ValueError):
            simple.trace_orbits(batch.copy(), integrator="does_not_exist", tmax=5e-5)


class TestSurfaceSampler:
    def test_surface_sampler_shape(self, sampler: simple.SurfaceSampler):
        batch = sampler.sample_surface_fieldline(8, s=0.4)
        assert batch.shape == (5, 8)
        assert batch.flags.c_contiguous


class TestParticleBatch:
    def test_initialize_from_samplers(self, vmec_file: str):
        batch = simple.ParticleBatch(12)
        batch.initialize_from_samplers(vmec_file, method="surface", s=0.35)

        assert batch.positions.shape == (5, 12)
        assert not np.allclose(batch.positions, 0.0)

        original = batch.positions[0, 0]
        batch.positions[0, 0] = original + 1.0
        assert batch.positions[0, 0] == original + 1.0

    def test_copy(self, vmec_file: str):
        batch = simple.ParticleBatch(6)
        batch.initialize_from_samplers(vmec_file, method="surface", s=0.5)

        clone = batch.copy()
        np.testing.assert_allclose(batch.positions, clone.positions)
        clone.positions[0, 0] += 0.1
        assert clone.positions[0, 0] != batch.positions[0, 0]


class TestTraceOrbits:
    def test_trace_returns_batch_results(self, vmec_file: str):
        results = run_trace_example(vmec_file, n_particles=4, tmax=1e-4)

        assert isinstance(results, simple.BatchResults)
        assert results.n_particles == 4
        assert results.final_positions.shape == (5, 4)
        assert results.loss_times.shape == (4,)
        assert results.tmax == pytest.approx(1e-4)

    def test_confined_and_lost_helpers(self, vmec_file: str):
        results = run_trace_example(vmec_file, n_particles=3, tmax=5e-5)

        confined = simple.get_confined(results)
        lost = simple.get_lost(results)

        assert confined.shape[0] == 5
        assert set(lost.keys()) == {"positions", "loss_times"}
        assert lost["positions"].shape[0] == 5


class TestParameterHelpers:
    def test_set_and_get_parameters(self, vmec_file: str):
        simple.load_vmec(vmec_file)

        original = simple.get_parameters("ntimstep")["ntimstep"]
        try:
            simple.set_parameters(ntimstep=2048)
            params = simple.get_parameters("ntimstep")
            assert params["ntimstep"] == 2048
        finally:
            simple.set_parameters(ntimstep=original)


class TestExampleModules:
    def test_fast_classification_example_disables_minkowski(self, vmec_file: str):
        result = classify_fast_example(vmec_file)

        assert result.j_parallel.shape == (16,)
        assert result.topology.shape == (16,)
        assert result.minkowski is None
