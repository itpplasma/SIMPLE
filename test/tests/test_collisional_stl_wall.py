#!/usr/bin/env python3
"""Seeded alpha slowing down and first STL hits through the production executable."""

from __future__ import annotations

import csv
import os
from pathlib import Path
import struct
import subprocess
import sys

from shapely import Polygon, points
from netCDF4 import Dataset
import numpy as np


NMARKERS = 192
TRACE_TIME = 5.0e-4
COIL_CURRENT_SCALE = 3.3
TRANSIT_TIME = 7.0e-6  # 2*pi*R_major/v0 at 35 keV, R_major about 1.4 m
RESULT_VARS = (
    "xstart_cart", "zstart", "zend", "times_lost", "orbit_exit_code", "wall_hit",
    "wall_hit_cart", "wall_hit_normal_cart", "wall_hit_angle_rad", "xend_cart",
)


def check_inside_lcfs(chartmap: Path, wout: Path) -> float:
    """Certify the launch ring on every chart slice against VMEC's LCFS."""
    with Dataset(wout) as ds:
        xm, xn = (np.asarray(ds[key][:]) for key in ("xm", "xn"))
        rmnc, zmns = (np.asarray(ds[key][-1]) for key in ("rmnc", "zmns"))
    with Dataset(chartmap) as ds:
        ir = int(np.argmin(abs(ds["rho"][:] - 0.25)))
        launch_rho = float(ds["rho"][ir])
        assert abs(launch_rho - 0.25) < 1e-6
        xyz = np.stack([np.asarray(ds[key][:, :, ir])
                        for key in ("x", "y", "z")], -1)
        xyz *= {"m": 1.0, "cm": 0.01}[ds["x"].units]
        zeta = np.asarray(ds["zeta"][:])
    theta = np.linspace(0, 2 * np.pi, 1024, endpoint=False)
    clearance = np.inf
    for phi, ring in zip(zeta, xyz):
        phase = xm[:, None] * theta - xn[:, None] * phi
        lcfs = Polygon(np.column_stack((rmnc @ np.cos(phase),
                                        zmns @ np.sin(phase))))
        launch = points(np.column_stack((np.hypot(ring[:, 0], ring[:, 1]),
                                         ring[:, 2])))
        assert lcfs.is_valid and np.all(lcfs.contains(launch)), "starts outside LCFS"
        clearance = min(clearance, np.min(lcfs.exterior.distance(launch)))
    assert clearance > 0.04, "launch ring lost its LCFS margin"
    print(f"Launch ring: minimum LCFS clearance {clearance:.6f} m")
    return launch_rho


def write_wall(chartmap: Path, output: Path) -> None:
    """Triangulate the boundary with a 5% inward margin inside the field chart."""
    with Dataset(chartmap) as ds:
        nfp = int(ds["num_field_periods"][...])
        surface = np.stack([np.asarray(ds[name][:, :, -1])
                            for name in ("x", "y", "z")], -1)
        axis = np.stack([np.asarray(ds[name][:, :, 0])
                         for name in ("x", "y", "z")], -1)
        units = ds["x"].units
    # Planar STL facets and the curved chart boundary do not coincide between
    # grid points. Leave a field-domain buffer so chart exits cannot mask hits.
    surface = axis + 0.95 * (surface - axis)
    scale = {"m": 1.0, "cm": 0.01}[units]
    # Chartmap angular grids exclude the periodic endpoint. Close both seams.
    surface = surface * scale
    rings = []
    for period in range(nfp):
        phi = period * 2 * np.pi / nfp
        rotation = np.array([[np.cos(phi), -np.sin(phi), 0],
                             [np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
        rings.append(surface @ rotation.T)
    vertices = np.concatenate(rings)
    nz, nt, _ = vertices.shape
    with output.open("wb") as stream:
        stream.write(b"NCSX chartmap outer boundary".ljust(80, b" "))
        stream.write(struct.pack("<I", 2 * nz * nt))
        for iz in range(nz):
            for it in range(nt):
                a = vertices[iz, it]
                b = vertices[(iz + 1) % nz, it]
                c = vertices[(iz + 1) % nz, (it + 1) % nt]
                d = vertices[iz, (it + 1) % nt]
                for triangle in ((a, b, c), (a, c, d)):
                    normal = np.cross(triangle[1] - triangle[0],
                                      triangle[2] - triangle[0])
                    normal /= np.linalg.norm(normal)
                    stream.write(struct.pack("<12fH", *normal,
                                             *np.asarray(triangle).ravel(), 0))


def run_case(base: Path, name: str, simple: Path, starts: np.ndarray,
             *, collisions: bool, density: float, npoiper2: int = 64) -> dict[str, np.ndarray]:
    work = base / name
    work.mkdir(exist_ok=True)
    (work / "results.nc").unlink(missing_ok=True)
    for filename in ("wall.stl", "chartmap.nc", "coils.simple", "wout.nc"):
        link = work / filename
        if not link.exists():
            link.symlink_to(base / filename)
    np.savetxt(work / "start.dat", starts, fmt="%.17e")
    (work / "simple.in").write_text(f"""&config
field_input = 'coils.simple'
coord_input = 'chartmap.nc'
netcdffile = 'wout.nc'
wall_input = 'wall.stl'
wall_units = 'm'
isw_field_type = 5
integmode = 0
startmode = 2
ntestpart = {len(starts)}
trace_time = {TRACE_TIME:.17e}
ntimstep = 101
npoiper2 = {npoiper2}
! Orbits use the coil field and chartmap; VMEC supplies equilibrium metadata.
multharm = 1
contr_pp = -1d10
relerr = 1d-9
facE_al = 1d2
deterministic = .true.
ran_seed = 541
swcoll = {'.true.' if collisions else '.false.'}
output_results_netcdf = .true.
/
&profiles
profile_type = 'power_series'
Te_scale = 1d3
Ti1_scale = 1d3
Ti2_scale = 1d3
ni1_scale = {density:.17e}
ni2_scale = {density:.17e}
/
""", encoding="ascii")
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE",
               GFORTRAN_UNBUFFERED_ALL="1")
    with (work / "run.log").open("w") as log:
        result = subprocess.run([str(simple)], cwd=work, env=env,
                                stdout=log, stderr=subprocess.STDOUT, timeout=360)
    assert result.returncode == 0, f"{name}: see {work / 'run.log'}"
    log = (work / "run.log").read_text()
    assert "Error in stost" not in log, f"{name}: collision error; see run.log"
    with Dataset(work / "results.nc") as ds:
        result = {key: np.asarray(ds[key][:]) for key in RESULT_VARS}
    assert not np.any(result["orbit_exit_code"] == 3), f"{name}: skipped markers"
    assert np.all(result["orbit_exit_code"] < 100), f"{name}: numerical exit"
    return result


def check_hits(result: dict[str, np.ndarray], name: str) -> np.ndarray:
    hit = result["wall_hit"] == 1
    assert np.count_nonzero(hit) >= 16, f"{name}: too few STL hits"
    assert np.all(result["orbit_exit_code"][hit] == 2), name
    assert np.all(result["times_lost"][hit] > 10 * TRANSIT_TIME), (
        f"{name}: prompt hit before ten toroidal transit times")
    assert np.all(result["times_lost"][hit] <= TRACE_TIME), name
    assert np.all(np.isfinite(result["zend"][hit])), name
    angle = result["wall_hit_angle_rad"][hit]
    assert np.all((angle >= 0) & (angle <= np.pi / 2)), name
    np.testing.assert_allclose(result["xend_cart"][hit],
                               result["wall_hit_cart"][hit], atol=1e-6, rtol=0)
    np.testing.assert_allclose(
        np.linalg.norm(result["wall_hit_normal_cart"][hit], axis=1), 1,
        atol=1e-12, rtol=0)
    print(f"{name}: {np.count_nonzero(hit)}/{NMARKERS} STL hits; "
          f"mean hit E/E0={np.mean(result['zend'][hit, 3] ** 2):.9f}; "
          f"hit times {np.min(result['times_lost'][hit]):.9e} to "
          f"{np.max(result['times_lost'][hit]):.9e} s")
    return hit


def main() -> None:
    simple, chartmap, coils, wout = (Path(arg).resolve() for arg in sys.argv[1:])
    base = Path.cwd() / "collisional_stl_wall"
    base.mkdir(exist_ok=True)
    for source, name in ((chartmap, "chartmap.nc"), (wout, "wout.nc")):
        link = base / name
        if not link.exists():
            link.symlink_to(source)
    # Scale only this test copy; keep the shared public coil fixture intact.
    data = np.loadtxt(coils, skiprows=1)
    data[:, 3] *= COIL_CURRENT_SCALE
    coil_copy = base / "coils.simple"
    coil_copy.unlink(missing_ok=True)
    np.savetxt(coil_copy, data, header=str(len(data)), comments="", fmt="%.17e")
    write_wall(chartmap, base / "wall.stl")
    # Explicit starts give every run the same ensemble independently of RNG use.
    launch_rho = check_inside_lcfs(chartmap, wout)
    with Dataset(chartmap) as ds:
        theta, zeta = (np.asarray(ds[key][:]) for key in ("theta", "zeta"))
    rng = np.random.default_rng(541)
    starts = np.column_stack((np.full(NMARKERS, launch_rho),
                              rng.choice(theta, NMARKERS),
                              rng.choice(zeta, NMARKERS),
                              np.ones(NMARKERS), rng.uniform(-0.2, 0.2, NMARKERS)))
    reference = run_case(base, "collisionless", simple, starts,
                         collisions=False, density=2.5e20)
    hit_ref = check_hits(reference, "collisionless")
    # Stored Fast-build reference; allow small platform-dependent orbit shifts.
    assert abs(np.count_nonzero(hit_ref) - 21) <= 2, "collisionless hit-count drift"
    np.testing.assert_allclose(np.mean(reference["times_lost"][hit_ref]),
                               2.3755234756525815e-4, rtol=0.03, atol=0)
    # Resolve the Cartesian segment interpolation independently of the stored
    # reference. Use known hits so this inexpensive subset covers real impacts.
    indices = np.flatnonzero(hit_ref)[:32]
    fine = run_case(base, "collisionless_fine", simple, starts[indices],
                    collisions=False, density=2.5e20, npoiper2=128)
    np.testing.assert_array_equal(fine["wall_hit"], 1)
    np.testing.assert_allclose(fine["times_lost"], reference["times_lost"][indices],
                               rtol=0.01, atol=1e-7)
    np.testing.assert_allclose(fine["wall_hit_cart"],
                               reference["wall_hit_cart"][indices],
                               rtol=0, atol=0.5)  # cm; below STL mesh resolution
    print("PASS: collisionless 64/128-interval segment convergence")
    dilute = run_case(base, "dilute", simple, starts,
                      collisions=True, density=0.0)
    # RK conserves the collisionless momentum; tolerance is well above relerr.
    np.testing.assert_allclose(reference["zend"][hit_ref, 3], 1,
                               rtol=0, atol=1e-7)
    # Exact vacuum traverses the collision path over the full wall trace.
    for key in ("wall_hit", "orbit_exit_code"):
        np.testing.assert_array_equal(dilute[key], reference[key], err_msg=key)
    for key in ("zend", "times_lost", "wall_hit_cart", "wall_hit_angle_rad"):
        np.testing.assert_allclose(dilute[key], reference[key],
                                   rtol=2e-6, atol=2e-8, err_msg=key)
    collisional = run_case(base, "collisional", simple, starts,
                           collisions=True, density=2.5e20)
    hit_coll = check_hits(collisional, "collisional")
    repeat = run_case(base, "repeat", simple, starts,
                      collisions=True, density=2.5e20)
    for key in RESULT_VARS:
        np.testing.assert_array_equal(repeat[key], collisional[key], err_msg=key)
    # Exclude prompt hits: these must have experienced several collision kicks.
    delayed = hit_coll & (collisional["times_lost"] > 20 * TRANSIT_TIME)
    assert np.count_nonzero(delayed) >= 16, "too few delayed collisional hits"
    energy = collisional["zend"][delayed, 3] ** 2
    assert np.mean(energy) < 0.99, "wall output lost the slowed-down state"
    assert np.count_nonzero(energy < 0.99) >= 16, "no resolved slowing down"
    with (base / "wall_hits.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("mode", "marker", "hit_energy_fraction", "hit_time_s",
                         "incidence_angle_rad"))
        for mode, result, hit in (("collisionless", reference, hit_ref),
                                  ("collisional", collisional, hit_coll)):
            for i in np.flatnonzero(hit):
                writer.writerow((mode, i + 1, result["zend"][i, 3] ** 2,
                                 result["times_lost"][i],
                                 result["wall_hit_angle_rad"][i]))
    print("PASS: collisionless limit, wall-hit slowing down, seeded repeat")
    print(f"CSV: {base / 'wall_hits.csv'}")


if __name__ == "__main__":
    main()
