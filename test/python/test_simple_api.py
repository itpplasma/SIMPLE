#!/usr/bin/env python3
"""
Tests for the cleaned SIMPLE Python API.
"""

from __future__ import annotations

import os
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


class TestConstants:
    """Test integration constants match Fortran."""

    def test_integration_constants(self):
        """Verify integrator constants match orbit_symplectic_base.f90"""
        assert pysimple.RK45 == 0
        assert pysimple.EXPL_IMPL_EULER == 1
        assert pysimple.IMPL_EXPL_EULER == 2
        assert pysimple.MIDPOINT == 3
        assert pysimple.GAUSS1 == 4
        assert pysimple.GAUSS2 == 5
        assert pysimple.GAUSS3 == 6
        assert pysimple.GAUSS4 == 7
        assert pysimple.LOBATTO3 == 15

    def test_default_values(self):
        """Verify default values match params.f90"""
        assert pysimple.DEFAULT_NS_S == 5
        assert pysimple.DEFAULT_NS_TP == 5
        assert pysimple.DEFAULT_MULTHARM == 5


class TestInitialization:
    """Test init() behavior."""

    def test_init_basic(self, vmec_file: str):
        """Basic initialization should succeed"""
        pysimple.init(vmec_file, deterministic=True)

    def test_init_with_parameters(self, vmec_file: str):
        """init() should accept parameter overrides"""
        pysimple.init(
            vmec_file,
            deterministic=True,
            trace_time=1e-3,
            ntestpart=100,
            npoiper2=64
        )

    def test_init_accepts_namelist_case(self, vmec_file: str):
        """init() should accept namelist-cased names such as facE_al."""
        pysimple.init(
            vmec_file,
            deterministic=True,
            trace_time=1e-4,
            ntestpart=1,
            facE_al=700.0,
            n_d=2,
            n_e=1,
        )
        assert float(pysimple.params.face_al) == 700.0

    def test_multiple_reinitializations(self, vmec_file: str):
        """Multiple init() calls should work without crash"""
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4)
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-3)
        pysimple.init(vmec_file, deterministic=True, trace_time=5e-5)

    def test_chartmap_init_selects_boozer_field(self):
        """Chart-map initialization must not attempt to load a VMEC wout."""
        chartmap_file = os.environ["SIMPLE_TEST_CHARTMAP"]
        pysimple.init(
            chartmap_file,
            deterministic=True,
            ntestpart=1,
            npoiper2=64,
            trace_time=1e-4,
        )
        assert int(pysimple._fortran_backend.velo_mod.isw_field_type) == 2
        result = pysimple.compute_canonical_frequencies(
            np.array([0.4, 0.7, 0.1, 1.0, 0.1]),
            n_periods=1,
            max_steps=10,
        )
        # A boundary stop is a physical loss (FREQ_ORBIT_LOST), never an
        # unexplained FREQ_INTEGRATOR_ERROR.
        assert result["status"] in {
            pysimple.FREQ_SUCCESS,
            pysimple.FREQ_ORBIT_LOST,
            pysimple.FREQ_MAX_STEPS,
        }


class TestIntegratorSelection:
    def test_alias_resolution(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True, trace_time=5e-5)
        particles = pysimple.sample_surface(2, s=0.42)

        pysimple.trace_parallel(particles.copy(), integrator="symplectic_midpoint")
        pysimple.trace_parallel(particles.copy(), integrator="rk45")

        with pytest.raises(ValueError):
            pysimple.trace_parallel(particles.copy(), integrator="does_not_exist")


class TestSampling:
    """Test particle sampling functions."""

    def test_surface_sampler_shape(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True)
        particles = pysimple.sample_surface(8, s=0.4)
        assert particles.shape == (5, 8)
        assert particles.flags.c_contiguous

    def test_sample_different_sizes(self, vmec_file: str):
        """Sample different numbers of particles"""
        pysimple.init(vmec_file, deterministic=True)
        for n in [1, 5, 20, 50]:
            particles = pysimple.sample_surface(n, s=0.5)
            assert particles.shape == (5, n)
            assert not np.allclose(particles, 0.0)

    def test_sample_volume(self, vmec_file: str):
        """Test volume sampling"""
        pysimple.init(vmec_file, deterministic=True)
        n_particles = 10
        particles = pysimple.sample_volume(n_particles, s_inner=0.2, s_outer=0.8)

        assert particles.shape == (5, n_particles)
        assert not np.allclose(particles, 0.0)
        # Check that s values are in expected range
        s_values = particles[0, :]
        assert np.all(s_values >= 0.2)
        assert np.all(s_values <= 0.8)


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
        # Explicitly set ntestpart to avoid interference
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=4)
        particles = pysimple.sample_surface(4, s=0.3)
        results = pysimple.trace_parallel(particles)

        assert isinstance(results, dict)
        assert results['final_positions'].shape == (5, 4)
        assert results['loss_times'].shape == (4,)
        assert results['trap_parameter'].shape == (4,)
        assert results['perpendicular_invariant'].shape == (4,)

    def test_trace_with_vmec_field(self, vmec_file: str):
        """Test tracing with VMEC field (isw_field_type=2)"""
        # Initialize with VMEC field type explicitly
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=3, isw_field_type=2)
        particles = pysimple.sample_surface(3, s=0.4)
        results = pysimple.trace_parallel(particles, integrator="midpoint")

        assert isinstance(results, dict)
        assert results['final_positions'].shape == (5, 3)
        assert results['loss_times'].shape == (3,)
        # Verify that trace completed (results are returned)
        assert 'trap_parameter' in results
        assert 'perpendicular_invariant' in results

    def test_trace_with_boozer_field(self, vmec_file: str):
        """Test tracing with Boozer canonical coordinates (isw_field_type=3)"""
        # Initialize with Boozer field type
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=2, isw_field_type=3)
        particles = pysimple.sample_surface(2, s=0.5)
        results = pysimple.trace_parallel(particles, integrator="midpoint")

        assert isinstance(results, dict)
        assert results['final_positions'].shape == (5, 2)
        assert results['loss_times'].shape == (2,)

    def test_confined_and_lost_filters(self, vmec_file: str):
        # Explicitly set ntestpart to avoid interference
        pysimple.init(vmec_file, deterministic=True, trace_time=5e-5, ntestpart=3)
        particles = pysimple.sample_surface(3, s=0.3)
        results = pysimple.trace_parallel(particles)

        confined_mask = results['loss_times'] >= 5e-5
        lost_mask = results['loss_times'] < 5e-5

        assert confined_mask.shape == (3,)
        assert lost_mask.shape == (3,)

    def test_different_integrators(self, vmec_file: str):
        """Test different integrator methods"""
        # Explicitly set ntestpart to avoid interference
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=2)
        particles = pysimple.sample_surface(2, s=0.5)

        # Test various integrators
        for integrator in ["midpoint", "rk45", "gauss2"]:
            results = pysimple.trace_parallel(particles.copy(), integrator=integrator)
            assert results['final_positions'].shape == (5, 2)
            assert results['loss_times'].shape == (2,)

    def test_trace_single_orbit(self, vmec_file: str):
        """Test trace_orbit for single particle"""
        # Re-init to avoid test interference
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=1)
        particles = pysimple.sample_surface(1, s=0.5)
        particle = particles[:, 0]

        result = pysimple.trace_orbit(particle, integrator="midpoint")
        assert result['final_position'].shape == (5,)
        assert isinstance(result['loss_time'], float)

    def test_trace_with_trajectory(self, vmec_file: str):
        """Test trace_orbit with trajectory output"""
        # Re-init to avoid test interference
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=1)
        particles = pysimple.sample_surface(1, s=0.5)
        particle = particles[:, 0]

        result = pysimple.trace_orbit(particle, integrator="midpoint", return_trajectory=True)
        assert 'trajectory' in result
        assert 'times' in result
        assert result['trajectory'].shape[0] == 5
        assert result['times'].shape[0] == result['trajectory'].shape[1]

    def test_canonical_frequencies(self, vmec_file: str):
        pysimple.init(
            vmec_file,
            deterministic=True,
            ntestpart=1,
            npoiper2=1024,
            trace_time=1e-3,
        )
        trapped = np.array([0.4, 0.7, 0.1, 1.0, 0.1])
        result = pysimple.compute_canonical_frequencies(
            trapped,
            integrator="midpoint",
            n_periods=2,
        )

        assert result["status"] == pysimple.FREQ_SUCCESS
        assert result["orbit_class"] == "trapped"
        assert result["n_periods"] == 2
        assert result["period"] > 0.0
        assert result["omega_b"] > 0.0
        assert result["omega_b"] == pytest.approx(2.0 * np.pi / result["period"])

    def test_canonical_frequency_validation(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True, ntestpart=1)
        particle = np.array([0.4, 0.7, 0.1, 1.0, 0.9])

        with pytest.raises(ValueError, match="n_periods"):
            pysimple.compute_canonical_frequencies(particle, n_periods=0)
        with pytest.raises(ValueError, match="max_steps"):
            pysimple.compute_canonical_frequencies(particle, max_steps=0)
        with pytest.raises(ValueError, match="shape"):
            pysimple.compute_canonical_frequencies(particle[:4])

    def test_trace_to_tip(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True, ntestpart=1, npoiper2=1024)
        particle = np.array([0.4, 0.7, 0.1, 1.0, 0.1])
        result = pysimple.trace_to_cut(particle, cut="tip", max_events=10)

        assert result["status"] == 0
        assert result["cut_type"] == "tip"
        assert result["state"].shape == (5,)
        assert abs(result["state"][4]) < 1e-6

    def test_trace_to_cut_validation(self, vmec_file: str):
        pysimple.init(vmec_file, deterministic=True, ntestpart=1)
        particle = np.array([0.4, 0.7, 0.1, 1.0, 0.9])
        with pytest.raises(ValueError, match="Unknown cut"):
            pysimple.trace_to_cut(particle, cut="poloidal")
        with pytest.raises(ValueError, match="symplectic"):
            pysimple.trace_to_cut(particle, integrator="rk45")

    def test_trace_parallel_preserves_skipped_passing_sentinel(self, vmec_file: str):
        """Deep-passing particles skipped by contr_pp must keep loss_time = -1."""
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=1)
        particle = np.array([[0.6], [0.0], [0.25 * np.pi], [1.0], [0.8]])

        trace_results = pysimple.trace_parallel(particle, integrator="midpoint")
        classify_results = pysimple.classify_parallel(particle, integrator="midpoint")

        assert trace_results["loss_times"][0] == pytest.approx(-1.0)
        assert classify_results["loss_times"][0] == pytest.approx(-1.0)
        np.testing.assert_allclose(trace_results["final_positions"][:, 0], particle[:, 0])
        np.testing.assert_allclose(classify_results["final_positions"][:, 0], particle[:, 0])

    def test_trace_orbit_uses_backend_loss_time_for_early_loss(self, vmec_file: str):
        """Single-particle tracing should report the backend loss time, not a padded zero."""
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-2, ntestpart=1)
        particle = np.array([0.6, 0.0, 0.25 * np.pi, 1.0, 0.2])

        single = pysimple.trace_orbit(
            particle,
            integrator="midpoint",
            return_trajectory=True,
        )

        pysimple.init(vmec_file, deterministic=True, trace_time=1e-2, ntestpart=1)
        batch = pysimple.trace_parallel(particle.reshape(5, 1), integrator="midpoint")

        assert single["loss_time"] == pytest.approx(batch["loss_times"][0])
        assert 0.0 < single["loss_time"] < pysimple.params.trace_time

        positive_times = np.flatnonzero(single["times"] > 0.0)
        assert positive_times.size > 0
        np.testing.assert_array_equal(
            np.isnan(single["times"][positive_times[-1] + 1 :]),
            np.ones(single["times"].shape[0] - positive_times[-1] - 1, dtype=bool),
        )
