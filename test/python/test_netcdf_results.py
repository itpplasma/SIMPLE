#!/usr/bin/env python3
"""
Tests for NetCDF results output (results.nc).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

try:
    import netCDF4 as nc
    HAS_NETCDF4 = True
except ImportError:
    HAS_NETCDF4 = False

python_dir = Path(__file__).parent.parent.parent / "python"
pytest.importorskip("pysimple", reason="pysimple module not available")

import sys
sys.path.insert(0, str(python_dir))

import pysimple


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
class TestNetCDFResultsOutput:
    """Test NetCDF results output functionality."""

    def test_results_netcdf_created(self, vmec_file: str, tmp_path: Path):
        """Test that results.nc is created when output_results_netcdf=True."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            pysimple.init(
                vmec_file,
                deterministic=True,
                trace_time=1e-5,
                ntestpart=8,
                output_results_netcdf=True,
                isw_field_type=3,
                contr_pp=-1e10,
                notrace_passing=0,
            )
            particles = pysimple.sample_surface(8, s=0.5)
            pysimple.trace_parallel(particles)
            pysimple.write_output()

            results_path = tmp_path / "results.nc"
            assert results_path.exists(), "results.nc was not created"

            with nc.Dataset(results_path, 'r') as ds:
                assert 'particle' in ds.dimensions
                assert 'xyz' in ds.dimensions
                assert 'phase' in ds.dimensions
                assert ds.dimensions['particle'].size == 8
                assert ds.dimensions['xyz'].size == 3
                assert ds.dimensions['phase'].size == 5

                assert 'times_lost' in ds.variables
                assert 'orbit_exit_code' in ds.variables
                assert 'boundary_event_radial_residual' in ds.variables
                assert 'boundary_event_time_width' in ds.variables
                assert 'zstart' in ds.variables
                assert 'zend' in ds.variables
                assert 'xstart_cart' in ds.variables
                assert 'xend_cart' in ds.variables
                assert 'trap_par' in ds.variables
                assert 'perp_inv' in ds.variables
                assert 'iclass' in ds.variables
                assert 'class_lost' in ds.variables
                assert np.all(np.isin(ds.variables['orbit_exit_code'][:], [0, 1, 2, 3]))
                exit_codes = ds.variables['orbit_exit_code'][:]
                event_residual = ds.variables['boundary_event_radial_residual'][:]
                event_width = ds.variables['boundary_event_time_width'][:]
                lcfs = exit_codes == 1
                assert np.all(event_residual[lcfs] >= 0.0)
                assert np.all(event_width[lcfs] >= 0.0)

                assert 'ntestpart' in ds.ncattrs()
                assert 'trace_time' in ds.ncattrs()
                assert ds.ntestpart == 8

        finally:
            os.chdir(original_dir)

    def test_no_results_when_disabled(self, vmec_file: str, tmp_path: Path):
        """Test that results.nc is NOT created when output_results_netcdf=False."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            pysimple.init(
                vmec_file,
                deterministic=True,
                trace_time=1e-5,
                ntestpart=4,
                output_results_netcdf=False,
                isw_field_type=3,
            )
            particles = pysimple.sample_surface(4, s=0.5)
            pysimple.trace_parallel(particles)
            pysimple.write_output()

            results_path = tmp_path / "results.nc"
            assert not results_path.exists(), \
                "results.nc should not be created when disabled"

        finally:
            os.chdir(original_dir)

    def test_cartesian_positions_at_known_surface(self, vmec_file: str, tmp_path: Path):
        """Test Cartesian positions are physically reasonable for s=0.5 surface."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            pysimple.init(
                vmec_file,
                deterministic=True,
                trace_time=1e-5,
                ntestpart=8,
                output_results_netcdf=True,
                isw_field_type=3,
            )
            particles = pysimple.sample_surface(8, s=0.5)
            pysimple.trace_parallel(particles)
            pysimple.write_output()

            with nc.Dataset(tmp_path / "results.nc", 'r') as ds:
                xstart = ds.variables['xstart_cart'][:]

                # For this VMEC file, R0 ~ 1022 cm
                # At s=0.5, particles should be within minor radius of axis
                # R = sqrt(x^2 + y^2) should be roughly 700-1300 cm
                r_cyl = np.sqrt(xstart[:, 0]**2 + xstart[:, 1]**2)
                assert np.all(r_cyl > 500), f"R too small: min={r_cyl.min():.1f}"
                assert np.all(r_cyl < 1500), f"R too large: max={r_cyl.max():.1f}"

                # Z should be within minor radius of midplane
                assert np.all(np.abs(xstart[:, 2]) < 500), \
                    f"Z out of range: max={np.abs(xstart[:, 2]).max():.1f}"

        finally:
            os.chdir(original_dir)

    def test_untraced_particles_have_zero_xend(
            self, vmec_file: str, tmp_path: Path):
        """Test that untraced particles (zend=0) have zero xend_cart."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            pysimple.init(
                vmec_file,
                deterministic=True,
                trace_time=1e-5,
                ntestpart=16,
                output_results_netcdf=True,
                isw_field_type=3,
            )
            particles = pysimple.sample_surface(16, s=0.5)
            pysimple.trace_parallel(particles)
            pysimple.write_output()

            with nc.Dataset(tmp_path / "results.nc", 'r') as ds:
                xend = ds.variables['xend_cart'][:]
                zend = ds.variables['zend'][:]

                for i in range(zend.shape[0]):
                    zend_is_zero = np.allclose(zend[i, :3], 0)
                    if zend_is_zero:
                        assert np.allclose(xend[i], 0), \
                            f"Particle {i}: untraced but xend_cart is nonzero"

        finally:
            os.chdir(original_dir)

    def test_wall_hits_use_authoritative_cartesian_hit(
            self, vmec_file: str, tmp_path: Path):
        """Test STL wall hits store xend_cart at wall_hit_cart."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            pysimple.init(
                vmec_file,
                deterministic=True,
                trace_time=1e-5,
                ntestpart=16,
                output_results_netcdf=True,
                isw_field_type=3,
            )
            particles = pysimple.sample_surface(16, s=0.5)
            pysimple.trace_parallel(particles)
            pysimple.write_output()

            with nc.Dataset(tmp_path / "results.nc", 'r') as ds:
                if 'wall_hit' not in ds.variables:
                    pytest.skip("wall-hit variables not available")

                wall_hit = ds.variables['wall_hit'][:]
                hit_indices = np.where(wall_hit == 1)[0]
                if hit_indices.size == 0:
                    pytest.skip("test run produced no STL wall hits")

                xend = ds.variables['xend_cart'][:]
                wall_hit_cart = ds.variables['wall_hit_cart'][:]
                assert np.allclose(xend[hit_indices], wall_hit_cart[hit_indices])

        finally:
            os.chdir(original_dir)
