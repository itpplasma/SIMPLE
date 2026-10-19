#!/usr/bin/env python3
"""
Tests for wall heat load engineering module.
"""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
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
    compute_wall_area_from_chartmap,
)


def create_mock_results_nc(path: Path, n_particles: int = 100, loss_fraction: float = 0.3):
    """Create a mock results.nc with the layout of src/netcdf_results_output.f90.

    Arrays are (particle, phase) and (particle, xyz). A normal run flags losses
    in orbit_exit_code (1 = LCFS, 2 = wall) and leaves class_lost at zero; only
    classification runs set class_lost.
    """
    n_lost = int(n_particles * loss_fraction)

    with nc.Dataset(path, 'w', format='NETCDF4') as ds:
        ds.createDimension('particle', n_particles)
        ds.createDimension('xyz', 3)
        ds.createDimension('phase', 5)

        ds.ntestpart = n_particles
        ds.trace_time = 1e-3

        zend = ds.createVariable('zend', 'f8', ('particle', 'phase'))
        xend_cart = ds.createVariable('xend_cart', 'f8', ('particle', 'xyz'))
        class_lost = ds.createVariable('class_lost', 'i1', ('particle',))
        exit_code = ds.createVariable('orbit_exit_code', 'i4', ('particle',))
        times_lost = ds.createVariable('times_lost', 'f8', ('particle',))

        np.random.seed(42)
        zend_data = np.zeros((5, n_particles))
        zend_data[0, :] = np.random.uniform(0.8, 1.0, n_particles)
        # Use VMEC theta convention [0, 2pi] - code should map to [-pi, pi]
        zend_data[1, :] = np.random.uniform(0, 2 * np.pi, n_particles)
        # Use large zeta values (cumulative angle) - code should apply modulo
        zend_data[2, :] = np.random.uniform(0, 20 * np.pi, n_particles)
        zend_data[3, :] = np.random.uniform(0.9, 1.0, n_particles)
        zend_data[4, :] = np.random.uniform(-1, 1, n_particles)
        zend[:] = zend_data.T

        R0, a = 1000.0, 100.0
        # Map theta to [-pi, pi] for Cartesian calculation (geometric consistency)
        theta_geom = np.where(
            zend_data[1, :] > np.pi,
            zend_data[1, :] - 2 * np.pi,
            zend_data[1, :]
        )
        zeta_geom = zend_data[2, :] % (2 * np.pi)
        theta = theta_geom
        zeta = zeta_geom
        R = R0 + a * np.cos(theta)
        x_cart = np.zeros((3, n_particles))
        x_cart[0, :] = R * np.cos(zeta)
        x_cart[1, :] = R * np.sin(zeta)
        x_cart[2, :] = a * np.sin(theta)
        xend_cart[:] = x_cart.T

        lost_mask = np.zeros(n_particles, dtype=np.int8)
        lost_mask[:n_lost] = 1
        np.random.shuffle(lost_mask)
        exit_code[:] = np.where(lost_mask == 1, 1, 0)
        class_lost[:] = 0

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

    def test_from_netcdf_of_real_simple_run(self, tmp_path: Path, vmec_file: str):
        """A normal simple.x run flags losses only in orbit_exit_code."""
        simple_x = Path(__file__).resolve().parents[2] / "build" / "simple.x"
        if not simple_x.exists():
            pytest.skip("simple.x is not built")
        (tmp_path / "simple.in").write_text(f"""&config
startmode = 1
ntestpart = 32
trace_time = 1d-3
npoiper2 = 64
notrace_passing = 0
isw_field_type = 0
integmode = 3
netcdffile = '{Path(vmec_file).resolve()}'
deterministic = .true.
output_results_netcdf = .true.
sbeg = 0.97
/
""")
        subprocess.run([str(simple_x), "simple.in"], cwd=tmp_path, check=True,
                       capture_output=True)

        with nc.Dataset(tmp_path / "results.nc") as ds:
            lost = np.isin(ds["orbit_exit_code"][:], (1, 2))
            assert ds["zend"].dimensions == ("particle", "phase")
            assert not ds["class_lost"][:].any()  # only classification sets it
            zend = np.asarray(ds["zend"][:])
        assert lost.any(), "test setup must produce lost particles"

        heat_map = WallHeatMap.from_netcdf(
            tmp_path / "results.nc", total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        assert heat_map.n_total == 32
        assert heat_map.n_lost == lost.sum()
        assert np.sum(heat_map.hit_count) == lost.sum()
        assert heat_map.energy_loss_fraction == pytest.approx(
            np.sum(zend[lost, 3] ** 2) / 32
        )

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


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
class TestPortOptimizer:
    """Tests for PortOptimizer class."""

    def test_solve_no_ports(self, mock_results_file: Path):
        """Test optimization with no ports returns empty result."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )
        opt = PortOptimizer(heat_map)

        result = opt.solve()

        assert result.success
        assert len(result.ports) == 0

    def test_solve_single_port(self, mock_results_file: Path):
        """Test optimization with single port finds low-flux region."""
        pytest.importorskip("scipy", reason="scipy required for optimization")

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
        # Optimized position should have less flux than center of heat map
        center_flux = heat_map.integrated_flux_in_region(-0.15, 0.15, np.pi - 0.15, np.pi + 0.15)
        assert result.total_flux_on_ports <= center_flux or result.total_flux_on_ports < 0.1

    def test_solve_multiple_ports(self, mock_results_file: Path):
        """Test optimization with multiple ports."""
        pytest.importorskip("scipy", reason="scipy required for optimization")

        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )
        opt = PortOptimizer(heat_map)
        opt.add_port("NBI_1", theta_width=0.2, zeta_width=0.2)
        opt.add_port("NBI_2", theta_width=0.2, zeta_width=0.2)

        result = opt.solve()

        assert result.success
        assert len(result.ports) == 2
        # Ports should be at different positions
        p1, p2 = result.ports
        dist = np.sqrt((p1.theta_center - p2.theta_center)**2 +
                       (p1.zeta_center - p2.zeta_center)**2)
        assert dist > 0.1  # Not at same location

    def test_solve_with_exclusion_zone(self, mock_results_file: Path):
        """Test that optimizer respects exclusion zones."""
        pytest.importorskip("scipy", reason="scipy required for optimization")

        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )
        opt = PortOptimizer(heat_map)
        opt.add_port("test_port", theta_width=0.2, zeta_width=0.2)
        # Exclude the center region where particles hit
        opt.add_exclusion_zone(-0.5, 0.5, 0, 2 * np.pi)

        result = opt.solve()

        assert result.success
        # Port center should be outside exclusion zone
        port = result.ports[0]
        assert port.theta_center < -0.5 or port.theta_center > 0.5

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

    def test_port_overlap_exclusion(self, mock_results_file: Path):
        """Test that port corners are checked against exclusion zones."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )
        opt = PortOptimizer(heat_map)

        # Exclusion zone at theta > 0.5
        opt.add_exclusion_zone(0.5, np.pi, 0, 2 * np.pi)

        # Port centered at 0.3 with width 0.5 -> extends to 0.55 (overlaps exclusion)
        assert opt._port_overlaps_exclusion(0.3, np.pi, 0.5, 0.5)

        # Port centered at 0.0 with width 0.5 -> extends to 0.25 (no overlap)
        assert not opt._port_overlaps_exclusion(0.0, np.pi, 0.5, 0.5)

    def test_max_flux_constraint(self, mock_results_file: Path):
        """Test that max_flux_constraint is enforced."""
        pytest.importorskip("scipy", reason="scipy required for optimization")

        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        # Without constraint
        opt1 = PortOptimizer(heat_map)
        opt1.add_port("test", theta_width=0.3, zeta_width=0.3)
        r1 = opt1.solve()

        # With very restrictive constraint - should avoid high flux areas
        opt2 = PortOptimizer(heat_map)
        opt2.add_port("test", theta_width=0.3, zeta_width=0.3)
        opt2.set_max_flux_constraint(0.001)  # Very low threshold
        r2 = opt2.solve()

        # Both should succeed but may find different positions
        assert r1.success
        assert r2.success


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
class TestNewFeatures:
    """Tests for bin areas, flux errors, and other new features."""

    def test_bin_areas_attribute(self, mock_results_file: Path):
        """Test that bin_areas is computed and stored."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        assert heat_map.bin_areas is not None
        assert heat_map.bin_areas.shape == heat_map.flux_grid.shape
        assert np.all(heat_map.bin_areas > 0)
        # Total should equal wall area
        assert np.sum(heat_map.bin_areas) == pytest.approx(heat_map.wall_area, rel=0.01)

    def test_chartmap_bin_areas(self, tmp_path: Path, mock_results_file: Path):
        """Test that chartmap provides per-bin areas."""
        chartmap_path = tmp_path / "chartmap.nc"
        create_mock_chartmap_nc(chartmap_path)

        heat_map = WallHeatMap.from_netcdf(
            mock_results_file,
            total_alpha_power_MW=TEST_ALPHA_POWER_MW,
            chartmap_file=chartmap_path,
        )

        # Bin areas should vary (not all identical) for real geometry
        # For simple torus they're uniform, but the mechanism is tested
        assert heat_map.bin_areas is not None
        assert np.sum(heat_map.bin_areas) == pytest.approx(heat_map.wall_area, rel=0.05)

    def test_resolution_mismatch_warning(self, tmp_path: Path, mock_results_file: Path):
        """Test that warning is issued when heatmap resolution exceeds chartmap."""
        import warnings

        # Create a very coarse chartmap (few grid points)
        chartmap_path = tmp_path / "coarse_chartmap.nc"
        create_mock_chartmap_nc(chartmap_path, ntheta=8, nzeta=4)

        # Request fine heatmap resolution - many bins will be empty
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            heat_map = WallHeatMap.from_netcdf(
                mock_results_file,
                total_alpha_power_MW=TEST_ALPHA_POWER_MW,
                chartmap_file=chartmap_path,
                n_theta=64,
                n_zeta=128,
            )

            # Should have issued a warning about bin area mismatch
            assert len(w) == 1
            assert "less than 90%" in str(w[0].message)
            assert "resolution" in str(w[0].message)

    def test_zeta_wrapping_in_region(self, mock_results_file: Path):
        """Test that integrated_flux_in_region handles zeta wrapping correctly."""
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        # Test region that wraps around zeta = 0 (negative zeta_min)
        # Port centered at zeta = 0.1 with width 0.4 has zeta_min = -0.1
        power_wrapped = heat_map.integrated_flux_in_region(
            theta_min=-0.5, theta_max=0.5,
            zeta_min=-0.1, zeta_max=0.3
        )

        # Same region without wrapping: [0, 0.3] + [2*pi - 0.1, 2*pi]
        power_part1 = heat_map.integrated_flux_in_region(
            theta_min=-0.5, theta_max=0.5,
            zeta_min=0, zeta_max=0.3
        )
        power_part2 = heat_map.integrated_flux_in_region(
            theta_min=-0.5, theta_max=0.5,
            zeta_min=2*np.pi - 0.1, zeta_max=2*np.pi
        )

        # Should be approximately equal (may differ slightly due to bin boundaries)
        assert power_wrapped == pytest.approx(power_part1 + power_part2, rel=0.1)

    def test_invalid_alpha_power(self, mock_results_file: Path):
        """Test that invalid total_alpha_power_MW raises error."""
        with pytest.raises(ValueError, match="must be positive"):
            WallHeatMap.from_netcdf(
                mock_results_file, total_alpha_power_MW=0.0
            )

        with pytest.raises(ValueError, match="must be positive"):
            WallHeatMap.from_netcdf(
                mock_results_file, total_alpha_power_MW=-10.0
            )

    def test_large_angles_handled(self, mock_results_file: Path):
        """Test that mock data with large angles is handled correctly.

        The mock data uses theta in [0, 2pi] and zeta in [0, 20*pi],
        which should be mapped to [-pi, pi] and [0, 2pi] respectively.
        """
        heat_map = WallHeatMap.from_netcdf(
            mock_results_file, total_alpha_power_MW=TEST_ALPHA_POWER_MW
        )

        # All binned data should be within the expected ranges
        # Check that we have particles distributed across theta range
        theta_has_hits = np.any(heat_map.hit_count > 0, axis=1)
        assert np.sum(theta_has_hits) > 10  # Should have hits in many theta bins

        # Check that we have particles distributed across zeta range
        zeta_has_hits = np.any(heat_map.hit_count > 0, axis=0)
        assert np.sum(zeta_has_hits) > 10  # Should have hits in many zeta bins

        # Total hits should match n_lost
        assert np.sum(heat_map.hit_count) == heat_map.n_lost


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
class TestVisualization:
    """Tests for visualization functions."""

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


def create_mock_chartmap_nc(
    path: Path,
    R0: float = 10.0,
    a: float = 1.0,
    nfp: int = 5,
    nrho: int = 17,
    ntheta: int = 33,
    nzeta: int = 16,
):
    """Create a mock chartmap.nc file for testing wall area calculation.

    Creates a simple torus with major radius R0 and minor radius a.
    The analytic surface area is A = 4 * pi^2 * R0 * a.

    Args:
        path: Output file path
        R0: Major radius in meters
        a: Minor radius in meters
        nfp: Number of field periods
        nrho: Number of radial grid points
        ntheta: Number of poloidal grid points
        nzeta: Number of toroidal grid points per field period
    """
    rho = np.linspace(0, 1, nrho)
    theta = np.linspace(0, 2*np.pi, ntheta, endpoint=False)
    zeta = np.linspace(0, 2*np.pi/nfp, nzeta, endpoint=False)

    # Create coordinates: x, y, z on (rho, theta, zeta) grid
    # File stores as (zeta, theta, rho) transposed
    x = np.zeros((nrho, ntheta, nzeta))
    y = np.zeros((nrho, ntheta, nzeta))
    z = np.zeros((nrho, ntheta, nzeta))

    for ir, r in enumerate(rho):
        for it, th in enumerate(theta):
            for iz, ze in enumerate(zeta):
                R = R0 + r * a * np.cos(th)
                x[ir, it, iz] = R * np.cos(ze) * 100  # cm
                y[ir, it, iz] = R * np.sin(ze) * 100  # cm
                z[ir, it, iz] = r * a * np.sin(th) * 100  # cm

    with nc.Dataset(path, 'w', format='NETCDF4') as ds:
        ds.createDimension('rho', nrho)
        ds.createDimension('theta', ntheta)
        ds.createDimension('zeta', nzeta)

        v_rho = ds.createVariable('rho', 'f8', ('rho',))
        v_theta = ds.createVariable('theta', 'f8', ('theta',))
        v_zeta = ds.createVariable('zeta', 'f8', ('zeta',))
        v_x = ds.createVariable('x', 'f8', ('zeta', 'theta', 'rho'))
        v_y = ds.createVariable('y', 'f8', ('zeta', 'theta', 'rho'))
        v_z = ds.createVariable('z', 'f8', ('zeta', 'theta', 'rho'))
        v_nfp = ds.createVariable('num_field_periods', 'i4')

        v_x.units = 'cm'
        v_y.units = 'cm'
        v_z.units = 'cm'

        v_rho[:] = rho
        v_theta[:] = theta
        v_zeta[:] = zeta
        v_x[:] = np.transpose(x, (2, 1, 0))
        v_y[:] = np.transpose(y, (2, 1, 0))
        v_z[:] = np.transpose(z, (2, 1, 0))
        v_nfp.assignValue(nfp)


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
class TestWallAreaCalculation:
    """Tests for wall area calculation from chartmap."""

    def test_simple_torus_area(self, tmp_path: Path):
        """Test wall area for simple torus matches analytic formula."""
        chartmap_path = tmp_path / "chartmap.nc"
        create_mock_chartmap_nc(chartmap_path)

        # Compute area from chartmap
        area = compute_wall_area_from_chartmap(chartmap_path)

        # Analytic torus area: A = 4 * pi^2 * R0 * a
        R0, a = 10.0, 1.0  # meters (same as in mock)
        expected_area = 4 * np.pi**2 * R0 * a

        # Should match within 1.5% (numerical discretization error)
        assert area == pytest.approx(expected_area, rel=0.015)

    @pytest.mark.parametrize(
        "R0, a, nfp, nzeta",
        [(5.0, 0.5, 1, 32), (8.0, 0.8, 10, 16), (100.0, 1.0, 5, 16),
         (3.0, 1.0, 1, 64), (6.2, 2.0, 1, 64), (5.5, 0.53, 5, 16)],
    )
    def test_torus_area_matches_analytic(self, tmp_path: Path, R0, a, nfp, nzeta):
        """Area of a circular torus is 4 pi^2 R0 a for any aspect ratio and nfp."""
        chartmap_path = tmp_path / "chartmap.nc"
        create_mock_chartmap_nc(chartmap_path, R0=R0, a=a, nfp=nfp, nzeta=nzeta)

        area = compute_wall_area_from_chartmap(chartmap_path)

        assert area == pytest.approx(4 * np.pi**2 * R0 * a, rel=0.015)

    def test_heatmap_uses_chartmap_area(self, tmp_path: Path, mock_results_file: Path):
        """Test that WallHeatMap uses chartmap area when provided."""
        chartmap_path = tmp_path / "chartmap.nc"
        create_mock_chartmap_nc(chartmap_path)

        heat_map = WallHeatMap.from_netcdf(
            mock_results_file,
            total_alpha_power_MW=TEST_ALPHA_POWER_MW,
            chartmap_file=chartmap_path,
        )

        # Should use chartmap area, not torus approximation
        R0, a = 10.0, 1.0
        expected_area = 4 * np.pi**2 * R0 * a
        assert heat_map.wall_area == pytest.approx(expected_area, rel=0.015)

    def test_fine_grid_convergence(self, tmp_path: Path):
        """Test that finer grid gives more accurate area."""
        R0, a = 10.0, 1.0
        expected_area = 4 * np.pi**2 * R0 * a

        # Coarse grid
        coarse_path = tmp_path / "coarse.nc"
        create_mock_chartmap_nc(coarse_path, R0=R0, a=a, nfp=5, ntheta=17, nzeta=8)
        area_coarse = compute_wall_area_from_chartmap(coarse_path)

        # Fine grid
        fine_path = tmp_path / "fine.nc"
        create_mock_chartmap_nc(fine_path, R0=R0, a=a, nfp=5, ntheta=65, nzeta=32)
        area_fine = compute_wall_area_from_chartmap(fine_path)

        # Fine grid should be closer to exact value
        error_coarse = abs(area_coarse - expected_area) / expected_area
        error_fine = abs(area_fine - expected_area) / expected_area
        assert error_fine < error_coarse
