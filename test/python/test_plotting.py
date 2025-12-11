#!/usr/bin/env python3
"""
Tests for pysimple.plotting module.
"""

from __future__ import annotations

from pathlib import Path
import tempfile

import numpy as np
import pytest

python_dir = Path(__file__).parent.parent.parent / "python"

import sys

sys.path.insert(0, str(python_dir))

from pysimple.plotting import (
    LossData,
    load_loss_data,
    load_slowing_down_curve,
    get_slowing_down_curve_path,
    compute_energy_confined_fraction,
    plot_energy_loss_vs_jperp,
    plot_confined_fraction,
)


class TestSlowingDownCurve:
    """Test slowing-down curve loading."""

    def test_get_bundled_path(self):
        """Test that bundled data file path is returned."""
        path = get_slowing_down_curve_path()
        # Should return a path-like object
        assert path is not None

    def test_load_bundled_curve(self):
        """Test loading bundled slowing-down curve."""
        times, energy = load_slowing_down_curve()
        assert len(times) == 100000
        assert len(energy) == 100000
        assert times[0] == 0.0
        assert np.isclose(times[-1], 1.0, rtol=1e-4)
        assert np.isclose(energy[0], 1.0, rtol=1e-4)
        assert energy[-1] < energy[0]  # Energy decreases

    def test_load_custom_path(self, tmp_path):
        """Test loading from custom path."""
        # Create a simple test file
        test_file = tmp_path / "test_curve.dat"
        test_data = np.column_stack([np.linspace(0, 1, 10), np.linspace(1, 0.1, 10)])
        np.savetxt(test_file, test_data)

        times, energy = load_slowing_down_curve(test_file)
        assert len(times) == 10
        assert len(energy) == 10


class TestLossData:
    """Test LossData loading and properties."""

    @pytest.fixture
    def sample_loss_data(self, tmp_path):
        """Create sample loss data files for testing."""
        n_particles = 100

        # times_lost.dat: index, t_loss, trap_par, s_start, J_perp, zend(1:5)
        times_lost = np.zeros((n_particles, 10))
        times_lost[:, 0] = np.arange(1, n_particles + 1)  # index
        times_lost[:, 1] = np.random.uniform(-1, 1, n_particles)  # t_loss
        times_lost[:, 2] = np.random.uniform(-1, 1, n_particles)  # trap_par
        times_lost[:, 3] = 0.25  # s_start
        times_lost[:, 4] = np.random.uniform(0, 1e-5, n_particles)  # J_perp
        times_lost[:, 8] = np.random.uniform(0.1, 1.0, n_particles)  # zend(4) = final_p

        # Make some particles "lost" (0 < t_loss < 1)
        lost_indices = np.random.choice(n_particles, 20, replace=False)
        times_lost[lost_indices, 1] = np.random.uniform(0.01, 0.9, 20)

        np.savetxt(tmp_path / "times_lost.dat", times_lost)

        # confined_fraction.dat: time, confpart_pass, confpart_trap, ntestpart
        n_times = 50
        conf_frac = np.zeros((n_times, 4))
        conf_frac[:, 0] = np.logspace(-5, 0, n_times)
        conf_frac[:, 1] = np.linspace(0.3, 0.25, n_times)  # passing
        conf_frac[:, 2] = np.linspace(0.6, 0.55, n_times)  # trapped
        conf_frac[:, 3] = n_particles

        np.savetxt(tmp_path / "confined_fraction.dat", conf_frac)

        # start.dat: s, theta, phi, v/v0, pitch
        start = np.zeros((n_particles, 5))
        start[:, 0] = 0.25
        start[:, 1] = np.random.uniform(0, 2 * np.pi, n_particles)
        start[:, 2] = np.random.uniform(0, 2 * np.pi, n_particles)
        start[:, 3] = 1.0
        start[:, 4] = np.random.uniform(-1, 1, n_particles)

        np.savetxt(tmp_path / "start.dat", start)

        return tmp_path

    def test_load_loss_data(self, sample_loss_data):
        """Test loading loss data from directory."""
        data = load_loss_data(sample_loss_data)

        assert isinstance(data, LossData)
        assert data.n_particles == 100
        assert len(data.loss_times) == 100
        assert len(data.trap_parameter) == 100
        assert len(data.perp_invariant) == 100
        assert len(data.final_p) == 100

    def test_loss_masks(self, sample_loss_data):
        """Test lost/confined/skipped masks."""
        data = load_loss_data(sample_loss_data)

        # Masks should be boolean arrays
        assert data.lost_mask.dtype == bool
        assert data.confined_mask.dtype == bool
        assert data.skipped_mask.dtype == bool

        # Masks should be mutually exclusive (mostly)
        # Note: particles with t_loss < 0 are skipped, t_loss >= trace_time are confined


class TestEnergyCalculations:
    """Test energy calculation functions."""

    @pytest.fixture
    def mock_loss_data(self):
        """Create mock LossData for testing."""
        n = 50
        return LossData(
            n_particles=n,
            loss_times=np.concatenate([np.full(10, -1), np.linspace(0.1, 0.9, 30), np.full(10, 1.0)]),
            trap_parameter=np.random.uniform(-1, 1, n),
            perp_invariant=np.random.uniform(0, 1e-5, n),
            start_s=np.full(n, 0.25),
            start_theta=np.random.uniform(0, 2 * np.pi, n),
            start_phi=np.random.uniform(0, 2 * np.pi, n),
            start_pitch=np.random.uniform(-1, 1, n),
            final_p=np.random.uniform(0.1, 1.0, n),
            time_grid=np.logspace(-5, 0, 100),
            confined_pass=np.linspace(0.3, 0.25, 100),
            confined_trap=np.linspace(0.6, 0.55, 100),
            trace_time=1.0,
        )

    def test_compute_energy_confined_fraction(self, mock_loss_data):
        """Test energy confined fraction calculation."""
        times, energy_frac = compute_energy_confined_fraction(mock_loss_data)

        assert len(times) == len(energy_frac)
        assert np.all(energy_frac >= 0)
        assert np.all(energy_frac <= 1)
        # Energy fraction should decrease over time as particles are lost
        assert energy_frac[0] >= energy_frac[-1]

    def test_compute_energy_with_slowing_down(self, mock_loss_data):
        """Test energy calculation with slowing-down curve."""
        sd_curve = load_slowing_down_curve()

        times, energy_frac = compute_energy_confined_fraction(
            mock_loss_data, slowing_down_curve=sd_curve
        )

        assert len(times) == len(energy_frac)
        assert np.all(energy_frac >= 0)
        assert np.all(energy_frac <= 1)


class TestPlotting:
    """Test plotting functions."""

    @pytest.fixture
    def mock_loss_data_pair(self):
        """Create mock LossData for COLL and NOCOLL runs."""
        n = 100

        def make_data(with_collisions):
            loss_times = np.concatenate([
                np.full(30, -1),  # skipped
                np.linspace(0.1, 0.9, 50),  # lost
                np.full(20, 1.0),  # confined
            ])
            if with_collisions:
                final_p = np.random.uniform(0.1, 0.8, n)
            else:
                final_p = np.ones(n)  # No collisions = full energy

            return LossData(
                n_particles=n,
                loss_times=loss_times,
                trap_parameter=np.random.uniform(-1, 1, n),
                perp_invariant=np.random.uniform(0, 1e-5, n),
                start_s=np.full(n, 0.25),
                start_theta=np.random.uniform(0, 2 * np.pi, n),
                start_phi=np.random.uniform(0, 2 * np.pi, n),
                start_pitch=np.random.uniform(-1, 1, n),
                final_p=final_p,
                time_grid=np.logspace(-5, 0, 100),
                confined_pass=np.linspace(0.3, 0.25, 100),
                confined_trap=np.linspace(0.6, 0.55, 100),
                trace_time=1.0,
            )

        return make_data(True), make_data(False)

    def test_plot_energy_loss_vs_jperp(self, mock_loss_data_pair, tmp_path):
        """Test energy loss vs J_perp plot generation."""
        pytest.importorskip("matplotlib")
        data_coll, data_nocoll = mock_loss_data_pair

        output_path = tmp_path / "energy_loss.png"
        fig = plot_energy_loss_vs_jperp(
            data_coll,
            data_nocoll,
            output_path=output_path,
            show=False,
        )

        assert output_path.exists()
        assert fig is not None

    def test_plot_energy_loss_with_slowing_down(self, mock_loss_data_pair, tmp_path):
        """Test energy loss plot with slowing-down curve."""
        pytest.importorskip("matplotlib")
        data_coll, data_nocoll = mock_loss_data_pair
        sd_curve = load_slowing_down_curve()

        output_path = tmp_path / "energy_loss_4curves.png"
        fig = plot_energy_loss_vs_jperp(
            data_coll,
            data_nocoll,
            slowing_down_curve=sd_curve,
            output_path=output_path,
            show=False,
        )

        assert output_path.exists()
        assert fig is not None

    def test_plot_confined_fraction(self, mock_loss_data_pair, tmp_path):
        """Test confined fraction plot generation."""
        pytest.importorskip("matplotlib")
        data_coll, _ = mock_loss_data_pair

        output_path = tmp_path / "confined_fraction.png"
        fig = plot_confined_fraction(
            data_coll,
            output_path=output_path,
            show=False,
        )

        assert output_path.exists()
        assert fig is not None

    def test_plot_confined_fraction_with_slowing_down(self, mock_loss_data_pair, tmp_path):
        """Test confined fraction plot with slowing-down curve."""
        pytest.importorskip("matplotlib")
        data_coll, _ = mock_loss_data_pair
        sd_curve = load_slowing_down_curve()

        output_path = tmp_path / "confined_fraction_theo.png"
        fig = plot_confined_fraction(
            data_coll,
            slowing_down_curve=sd_curve,
            output_path=output_path,
            show=False,
        )

        assert output_path.exists()
        assert fig is not None
