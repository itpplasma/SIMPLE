#!/usr/bin/env python3
"""
Tests for wall heat load engineering module.
"""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

import numpy as np
import pytest

try:
    import netCDF4 as nc
    HAS_NETCDF4 = True
except ImportError:
    HAS_NETCDF4 = False

python_dir = Path(__file__).parent.parent.parent / "python"

import sys
sys.path.insert(0, str(python_dir))

from pysimple.engineering import (
    WallHeatMap,
    Port,
    PortOptimizer,
    OptimizationResult,
    plot_heat_flux_2d,
)


def create_mock_results_nc(path: Path, n_particles: int = 100, loss_fraction: float = 0.3):
    """Create a mock results.nc file for testing."""
    n_lost = int(n_particles * loss_fraction)

    with nc.Dataset(path, 'w', format='NETCDF4') as ds:
        ds.createDimension('particle', n_particles)
        ds.createDimension('xyz', 3)
        ds.createDimension('phase', 5)

        ds.ntestpart = n_particles
        ds.trace_time = 1e-3

        zend = ds.createVariable('zend', 'f8', ('phase', 'particle'))
        xend_cart = ds.createVariable('xend_cart', 'f8', ('xyz', 'particle'))
        class_lost = ds.createVariable('class_lost', 'i1', ('particle',))
        times_lost = ds.createVariable('times_lost', 'f8', ('particle',))

        np.random.seed(42)
        zend_data = np.zeros((5, n_particles))
        zend_data[0, :] = np.random.uniform(0.8, 1.0, n_particles)
        zend_data[1, :] = np.random.uniform(-np.pi, np.pi, n_particles)
        zend_data[2, :] = np.random.uniform(0, 2 * np.pi, n_particles)
        zend_data[3, :] = np.random.uniform(0.9, 1.0, n_particles)
        zend_data[4, :] = np.random.uniform(-1, 1, n_particles)
        zend[:] = zend_data

        R0, a = 1000.0, 100.0
        theta = zend_data[1, :]
        zeta = zend_data[2, :]
        R = R0 + a * np.cos(theta)
        x_cart = np.zeros((3, n_particles))
        x_cart[0, :] = R * np.cos(zeta)
        x_cart[1, :] = R * np.sin(zeta)
        x_cart[2, :] = a * np.sin(theta)
        xend_cart[:] = x_cart

        lost_mask = np.zeros(n_particles, dtype=np.int8)
        lost_mask[:n_lost] = 1
        np.random.shuffle(lost_mask)
        class_lost[:] = lost_mask

        t_lost = np.full(n_particles, -1.0)
        t_lost[lost_mask == 1] = np.random.uniform(0, 1e-3, n_lost)
        times_lost[:] = t_lost


@pytest.fixture
def mock_results_file(tmp_path: Path) -> Path:
    """Create a temporary mock results.nc file."""
    results_path = tmp_path / "results.nc"
    create_mock_results_nc(results_path, n_particles=100, loss_fraction=0.3)
    return results_path


# Default alpha power for tests: 600 MW (typical for 3 GW fusion reactor)
TEST_ALPHA_POWER_MW = 600.0


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
class TestWallHeatMap:
    """Tests for WallHeatMap class."""

    def test_load_from_netcdf(self, mock_results_file: Path):
        """Test loading heat map from NetCDF results."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        assert heat_map.n_total == 100
        assert heat_map.n_lost == 30
        assert heat_map.loss_fraction == pytest.approx(0.3, rel=0.1)
        assert heat_map.trace_time == pytest.approx(1e-3)
        assert heat_map.total_alpha_power == TEST_ALPHA_POWER_MW

    def test_flux_grid_shape(self, mock_results_file: Path):
        """Test flux grid has correct shape."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file,
            total_alpha_power_MW=TEST_ALPHA_POWER_MW,
            n_theta=32,
            n_zeta=64,
        )

        assert heat_map.flux_grid.shape == (32, 64)
        assert heat_map.hit_count.shape == (32, 64)
        assert len(heat_map.theta_edges) == 33
        assert len(heat_map.zeta_edges) == 65

    def test_lost_power_energy_weighted(self, mock_results_file: Path):
        """Test lost power uses energy weighting, not just particle count."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        # Energy loss fraction accounts for slowing down (p < 1)
        # In mock data, p = v/v0 is in [0.9, 1.0], so energy = p^2 < 1
        assert heat_map.energy_loss_fraction <= heat_map.loss_fraction
        assert heat_map.energy_loss_fraction > 0

        # Lost power = energy_loss_fraction * total_alpha_power
        expected_lost = heat_map.energy_loss_fraction * TEST_ALPHA_POWER_MW
        assert heat_map.lost_power == pytest.approx(expected_lost, rel=1e-10)
        assert heat_map.peak_flux >= 0
        assert heat_map.mean_flux >= 0

    def test_energy_grid_sums_correctly(self, mock_results_file: Path):
        """Test energy_grid sums to total energy lost."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        # Sum of energy_grid = total energy weight = n_total * energy_loss_fraction
        total_energy = np.sum(heat_map.energy_grid)
        expected = heat_map.n_total * heat_map.energy_loss_fraction
        assert total_energy == pytest.approx(expected, rel=1e-10)

    def test_hit_count_matches_lost(self, mock_results_file: Path):
        """Test total hit count equals number of lost particles."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        assert np.sum(heat_map.hit_count) == heat_map.n_lost

    def test_flux_at_location(self, mock_results_file: Path):
        """Test flux lookup at specific location."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        flux = heat_map.flux_at(0.0, np.pi)
        assert isinstance(flux, float)
        assert flux >= 0

    def test_integrated_flux_in_region(self, mock_results_file: Path):
        """Test integrated flux in a region."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        # Full region should capture all lost power
        full_power = heat_map.integrated_flux_in_region(
            -np.pi, np.pi, 0, 2 * np.pi
        )
        assert full_power == pytest.approx(heat_map.lost_power, rel=0.01)

        # Half region should be less
        half_power = heat_map.integrated_flux_in_region(
            -np.pi, np.pi, 0, np.pi
        )
        assert half_power < heat_map.lost_power

    def test_summary_string(self, mock_results_file: Path):
        """Test summary string generation."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        summary = heat_map.summary()
        assert "30/100" in summary
        assert "MW" in summary
        assert "Alpha power" in summary
        assert "Particle loss fraction" in summary
        assert "Energy loss fraction" in summary

    def test_no_lost_particles(self, tmp_path: Path):
        """Test handling of results with no lost particles."""
        results_path = tmp_path / "no_lost.nc"
        create_mock_results_nc(results_path, n_particles=50, loss_fraction=0.0)

        heat_map = WallHeatMap.from_netcdf(
            results_path, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        assert heat_map.n_lost == 0
        assert heat_map.lost_power == 0
        assert heat_map.energy_loss_fraction == 0
        assert heat_map.peak_flux == 0
        assert np.all(heat_map.flux_grid == 0)
        assert np.all(heat_map.energy_grid == 0)


class TestPort:
    """Tests for Port dataclass."""

    def test_port_bounds(self):
        """Test port boundary calculations."""
        port = Port(
            name="test",
            theta_center=0.5,
            zeta_center=1.0,
            theta_width=0.2,
            zeta_width=0.4,
        )

        assert port.theta_min == pytest.approx(0.4)
        assert port.theta_max == pytest.approx(0.6)
        assert port.zeta_min == pytest.approx(0.8)
        assert port.zeta_max == pytest.approx(1.2)


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
class TestPortOptimizer:
    """Tests for PortOptimizer class."""

    def test_add_port(self, mock_results_file: Path):
        """Test adding ports to optimizer."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )
        opt = PortOptimizer(heat_map)

        opt.add_port("NBI_1", theta_width=0.3, zeta_width=0.2)
        opt.add_port("diag_1", theta_width=0.1, zeta_width=0.1)

        assert len(opt._ports) == 2

    def test_solve_no_ports(self, mock_results_file: Path):
        """Test optimization with no ports returns empty result."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )
        opt = PortOptimizer(heat_map)

        result = opt.solve()

        assert result.success
        assert len(result.ports) == 0

    @pytest.mark.skipif(
        True,  # Skip by default - requires scipy
        reason="Optimization test requires scipy"
    )
    def test_solve_single_port(self, mock_results_file: Path):
        """Test optimization with single port."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )
        opt = PortOptimizer(heat_map)
        opt.add_port("test_port", theta_width=0.3, zeta_width=0.3)

        result = opt.solve()

        assert result.success
        assert len(result.ports) == 1
        assert result.ports[0].name == "test_port"
        assert -np.pi <= result.ports[0].theta_center <= np.pi
        assert 0 <= result.ports[0].zeta_center <= 2 * np.pi

    def test_exclusion_zone(self, mock_results_file: Path):
        """Test adding exclusion zones."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )
        opt = PortOptimizer(heat_map)

        opt.add_exclusion_zone(-0.5, 0.5, 0, np.pi)

        assert len(opt._exclusion_zones) == 1
        assert opt._in_exclusion_zone(0.0, 0.5)
        assert not opt._in_exclusion_zone(1.0, 0.5)


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
class TestVisualization:
    """Tests for visualization functions."""

    def test_plot_heat_flux_2d_no_display(self, mock_results_file: Path):
        """Test 2D plot creation without display."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        fig, ax = plt.subplots()
        result_ax = plot_heat_flux_2d(heat_map, ax=ax, show_colorbar=True)

        assert result_ax is ax
        plt.close(fig)

    def test_plot_with_ports(self, mock_results_file: Path):
        """Test 2D plot with port overlays."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )
        ports = [
            Port("NBI", 0.0, np.pi, 0.3, 0.2),
            Port("diag", 1.0, 2.0, 0.1, 0.1),
        ]

        fig, ax = plt.subplots()
        plot_heat_flux_2d(heat_map, ax=ax, ports=ports)

        assert len(ax.patches) == 2
        plt.close(fig)


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
class TestIntegration:
    """Integration tests with real SIMPLE output."""

    @pytest.fixture
    def vmec_file(self) -> str:
        """Get path to test VMEC file."""
        vmec_path = Path(__file__).parent.parent / "data" / "wout_vmec_simple.nc"
        if not vmec_path.exists():
            pytest.skip("Test VMEC file not available")
        return str(vmec_path)

    def test_full_workflow_mock(self, mock_results_file: Path):
        """Test full workflow with mock data."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        assert heat_map.n_lost > 0
        assert heat_map.lost_power > 0
        assert heat_map.energy_loss_fraction > 0

        opt = PortOptimizer(heat_map)
        opt.add_port("test", theta_width=0.2, zeta_width=0.2)

        summary = heat_map.summary()
        assert "lost" in summary.lower()
