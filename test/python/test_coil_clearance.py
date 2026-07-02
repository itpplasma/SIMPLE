#!/usr/bin/env python3
"""Tests for coil-clearance geometry helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

try:
    import netCDF4 as nc
    HAS_NETCDF4 = True
except ImportError:
    HAS_NETCDF4 = False

from pysimple.coil_clearance import (
    CoilClearanceConstraint,
    CoilSegments,
    WallSurface,
    clearance_map_on_heatmap_grid,
    load_coil_points,
    min_distance_points_to_segments,
    split_coil_polylines,
)


def write_coil_file(path: Path) -> None:
    data = np.asarray(
        [
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 1.0],
            [0.0, 2.0, 0.0, 1.0],
        ]
    )
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{data.shape[0]}\n")
        np.savetxt(f, data)


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
def write_chartmap(path: Path, *, nfp: int = 2) -> None:
    theta = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi / nfp, 6, endpoint=False)
    rho = np.asarray([0.0, 1.0])

    with nc.Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("zeta", zeta.size)
        ds.createDimension("theta", theta.size)
        ds.createDimension("rho", rho.size)

        ds.createVariable("theta", "f8", ("theta",))[:] = theta
        ds.createVariable("zeta", "f8", ("zeta",))[:] = zeta
        ds.createVariable("rho", "f8", ("rho",))[:] = rho
        ds.createVariable("num_field_periods", "i4").assignValue(nfp)

        x = ds.createVariable("x", "f8", ("zeta", "theta", "rho"))
        y = ds.createVariable("y", "f8", ("zeta", "theta", "rho"))
        z = ds.createVariable("z", "f8", ("zeta", "theta", "rho"))

        for iz, ze in enumerate(zeta):
            for it, th in enumerate(theta):
                radius_cm = 100.0 + 10.0 * np.cos(th)
                x[iz, it, 1] = radius_cm * np.cos(ze)
                y[iz, it, 1] = radius_cm * np.sin(ze)
                z[iz, it, 1] = 10.0 * np.sin(th)


def test_load_and_split_coil_file(tmp_path: Path):
    coil_file = tmp_path / "coils.dat"
    write_coil_file(coil_file)

    points, currents = load_coil_points(coil_file)
    polylines = split_coil_polylines(points, currents)
    coils = CoilSegments.from_file(coil_file)

    assert points.shape == (5, 3)
    assert len(polylines) == 2
    assert coils.p0.shape == (2, 3)
    assert coils.p1.shape == (2, 3)


def test_min_distance_points_to_segments():
    points = np.asarray([[0.5, 1.0, 0.0], [2.0, 0.0, 0.0]])
    p0 = np.asarray([[0.0, 0.0, 0.0]])
    p1 = np.asarray([[1.0, 0.0, 0.0]])

    distances = min_distance_points_to_segments(points, p0, p1)

    assert distances[0] == pytest.approx(1.0)
    assert distances[1] == pytest.approx(1.0)


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
def test_wall_surface_period_rotation(tmp_path: Path):
    chartmap = tmp_path / "chartmap.nc"
    write_chartmap(chartmap, nfp=2)
    wall = WallSurface.from_chartmap(chartmap)

    xyz_period0 = wall.xyz(np.asarray([0.0]), np.asarray([0.0]))[0]
    xyz_period1 = wall.xyz(np.asarray([0.0]), np.asarray([np.pi]))[0]

    assert xyz_period0[0] == pytest.approx(1.1)
    assert xyz_period0[1] == pytest.approx(0.0)
    assert xyz_period1[0] == pytest.approx(-1.1)
    assert xyz_period1[1] == pytest.approx(0.0, abs=1e-14)


@pytest.mark.skipif(not HAS_NETCDF4, reason="netCDF4 not available")
def test_coil_clearance_constraint_and_map(tmp_path: Path):
    chartmap = tmp_path / "chartmap.nc"
    write_chartmap(chartmap, nfp=1)
    wall = WallSurface.from_chartmap(chartmap)

    coils = CoilSegments(
        p0=np.asarray([[1.15, 0.0, -0.5]]),
        p1=np.asarray([[1.15, 0.0, 0.5]]),
    )
    constraint = CoilClearanceConstraint(
        wall=wall,
        coils=coils,
        min_clearance_m=0.1,
    )

    assert constraint.port_violates(0.0, 0.0, 0.1, 0.1)
    assert not constraint.port_violates(0.0, np.pi, 0.1, 0.1)

    clearance = clearance_map_on_heatmap_grid(
        wall,
        coils,
        np.asarray([0.0]),
        np.asarray([0.0, np.pi]),
    )
    assert clearance.shape == (1, 2)
    assert clearance[0, 0] < clearance[0, 1]
