#!/usr/bin/env python3
"""Reduced ensembles; the unchanged 192-marker test remains the full-tier oracle."""

from math import comb
from pathlib import Path
import os
import subprocess
import sys
from time import perf_counter

from netCDF4 import Dataset
import numpy as np

import test_collisional_stl_wall as full

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'python'))
from pysimple.wall_heat_load import TriangleHeatLoad


def check_hits(result, name):
    hit = result["wall_hit"] == 1
    assert np.any(hit), f"{name}: no STL hits"
    assert np.all(result["orbit_exit_code"][hit] == 2), name
    assert np.all(result["times_lost"][hit] > 10 * full.TRANSIT_TIME), name
    assert np.all(result["times_lost"][hit] <= full.TRACE_TIME), name
    assert np.all(np.isfinite(result["zend"][hit])), name
    angle = result["wall_hit_angle_rad"][hit]
    assert np.all((angle >= 0) & (angle <= np.pi / 2)), name
    np.testing.assert_allclose(
        result["xend_cart"][hit], result["wall_hit_cart"][hit], atol=1e-6, rtol=0
    )
    np.testing.assert_allclose(
        np.linalg.norm(result["wall_hit_normal_cart"][hit], axis=1),
        1,
        atol=1e-12,
        rtol=0,
    )
    print(f"{name}: {hit.sum()}/{len(hit)} STL hits", flush=True)
    return hit


def write_case(base, case):
    name, ensemble, collisions, intervals = case
    work = base / name
    work.mkdir(exist_ok=True)
    (work / "results.nc").unlink(missing_ok=True)
    for filename in ("wall.stl", "chartmap.nc", "coils.simple", "wout.nc"):
        link = work / filename
        if not link.exists():
            link.symlink_to(base / filename)
    np.savetxt(work / "start.dat", ensemble, fmt="%.17e")
    density = 0.0 if name == "dilute" else 2.5e20
    (work / "simple.in").write_text(
        f"""&config
field_input = 'coils.simple'
coord_input = 'chartmap.nc'
netcdffile = 'wout.nc'
wall_input = 'wall.stl'
wall_units = 'm'
isw_field_type = 5
integmode = 0
startmode = 2
ntestpart = {len(ensemble)}
trace_time = {full.TRACE_TIME:.17e}
ntimstep = 101
npoiper2 = {intervals}
multharm = 1
contr_pp = -1d10
relerr = 1d-9
facE_al = 1d2
deterministic = .true.
ran_seed = 541
swcoll = {".true." if collisions else ".false."}
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
""",
        encoding="ascii",
    )
    return work


def main():
    begin = perf_counter()
    simple, chartmap, coils, wout = (Path(arg).resolve() for arg in sys.argv[1:])
    base = Path.cwd() / "collisional_stl_wall_fast"
    base.mkdir(exist_ok=True)
    for source, name in ((chartmap, "chartmap.nc"), (wout, "wout.nc")):
        link = base / name
        if not link.exists():
            link.symlink_to(source)
    data = np.loadtxt(coils, skiprows=1)
    data[:, 3] *= full.COIL_CURRENT_SCALE
    np.savetxt(
        base / "coils.simple", data, header=str(len(data)), comments="", fmt="%.17e"
    )
    full.write_wall(chartmap, base / "wall.stl")
    rho = full.check_inside_lcfs(chartmap, wout)
    with Dataset(chartmap) as ds:
        theta, zeta = (np.asarray(ds[key][:]) for key in ("theta", "zeta"))
    rng = np.random.default_rng(541)
    # Subsample the original ensemble, without changing the physical duration,
    # fields, integration tolerances, collision density or wall resolution.
    starts = np.column_stack(
        (
            np.full(full.NMARKERS, rho),
            rng.choice(theta, full.NMARKERS),
            rng.choice(zeta, full.NMARKERS),
            np.ones(full.NMARKERS),
            rng.uniform(-0.2, 0.2, full.NMARKERS),
        )
    )
    ref_indices = np.linspace(0, full.NMARKERS - 1, 32, dtype=int)
    coll_indices = np.linspace(0, full.NMARKERS - 1, 64, dtype=int)
    # Stored collisionless impacts in the 32-marker subset; replay all five
    # at doubled segment resolution.
    fine_indices = np.array([43, 98, 117, 166, 172])
    cases = (
        ("collisionless", starts[ref_indices], False, 64),
        ("dilute", starts[ref_indices], True, 64),
        ("collisional", starts[coll_indices], True, 64),
        ("repeat", starts[coll_indices], True, 64),
        ("collisionless_fine", starts[fine_indices], False, 128),
    )

    workdirs = [write_case(base, case) for case in cases]
    (base / "simple.in").write_text((workdirs[0] / "simple.in").read_text())
    # One production-library driver initializes the field once. OpenMP uses
    # four threads for both field construction and per-marker integration.
    env = dict(
        os.environ,
        OMP_NUM_THREADS="4",
        OMP_DYNAMIC="FALSE",
        GFORTRAN_UNBUFFERED_ALL="1",
    )
    with (base / "run.log").open("w") as log:
        process = subprocess.run(
            [str(simple), *map(str, workdirs)],
            cwd=base,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=120,
        )
    assert process.returncode == 0, f"see {base / 'run.log'}"
    assert "Error in stost" not in (base / "run.log").read_text()
    results = []
    for work in workdirs:
        with Dataset(work / "results.nc") as ds:
            result = {key: np.asarray(ds[key][:]) for key in full.RESULT_VARS}
            energy = ds['kinetic_energy_eV']
            assert energy.units == 'eV'
            np.testing.assert_allclose(
                energy[:], ds.particle_energy_eV * result['zend'][:, 3] ** 2,
                rtol=2e-15, atol=0)
            assert energy.filters()['zlib']
            np.testing.assert_array_equal(ds['marker_weight'][:], 1)
            assert ds['marker_weight'].source.startswith('equal_marker_default;')
            assert 'chord' in ds.wall_incidence_representation
            assert ds['wall_hit_triangle'].start_index == 0
            assert 'NI1_SCALE' in ds.profile_configuration.upper()
            np.testing.assert_array_equal(ds.background_ion_charge_number, [1, 1])

        with Dataset(work / 'weighted-results.nc') as weighted:
            weights = np.ones(len(result['wall_hit']))
            weights[0] = 2
            np.testing.assert_array_equal(weighted['marker_weight'][:], weights)
            assert weighted['marker_weight'].source == 'supplied_relative_weights'
        hit = result['wall_hit'] == 1
        for filename, source_weights in (
                ('results.nc', np.ones(len(hit))),
                ('weighted-results.nc', weights)):
            loads = TriangleHeatLoad.from_netcdf(work / filename, work / 'wall.stl', 1)
            expected = np.sum(source_weights[hit] * result['zend'][hit, 3] ** 2)
            expected *= 1e6 / np.sum(source_weights)
            np.testing.assert_allclose(loads.deposited_power_W, expected, rtol=2e-15)
        assert not np.any(result["orbit_exit_code"] == 3), str(work)
        assert np.all(result["orbit_exit_code"] < 100), str(work)
        results.append(result)
    reference, dilute, collisional, repeat, fine = results
    hit_ref = check_hits(reference, "collisionless")
    hit_coll = check_hits(collisional, "collisional")
    # Numerical drift budgets are unchanged; the stored mean is recomputed for
    # this exact deterministic subset, rather than widened for a smaller sample.
    assert abs(np.count_nonzero(hit_ref) - 5) <= 2, "collisionless hit-count drift"
    np.testing.assert_allclose(
        np.mean(reference["times_lost"][hit_ref]),
        2.7560491342032043e-4,
        rtol=0.03,
        atol=0,
    )
    matched = np.searchsorted(ref_indices, fine_indices)
    np.testing.assert_array_equal(reference["wall_hit"][matched], 1)
    np.testing.assert_array_equal(fine["wall_hit"], 1)
    np.testing.assert_allclose(
        fine["times_lost"], reference["times_lost"][matched], rtol=0.01, atol=1e-7
    )
    np.testing.assert_allclose(
        fine["wall_hit_cart"], reference["wall_hit_cart"][matched], rtol=0, atol=0.5
    )
    np.testing.assert_allclose(reference["zend"][hit_ref, 3], 1, rtol=0, atol=1e-7)
    for key in ("wall_hit", "orbit_exit_code"):
        np.testing.assert_array_equal(dilute[key], reference[key], err_msg=key)
    for key in ("zend", "times_lost", "wall_hit_cart", "wall_hit_angle_rad"):
        np.testing.assert_allclose(
            dilute[key], reference[key], rtol=2e-6, atol=2e-8, err_msg=key
        )
    for key in full.RESULT_VARS:
        np.testing.assert_array_equal(repeat[key], collisional[key], err_msg=key)
    delayed = hit_coll & (collisional["times_lost"] > 20 * full.TRANSIT_TIME)
    energy = collisional["zend"][delayed, 3] ** 2
    # Exact one-sided sign test of median E/E0 >= 0.99: under the null, each
    # independent delayed impact has probability <= 1/2 of being below 0.99.
    # P[Binomial(n, 1/2) >= k] <= 0.05 requires at least five slowed impacts;
    # it automatically tightens if any delayed impacts fail to slow down.
    n = len(energy)
    k = np.count_nonzero(energy < 0.99)
    p_value = sum(comb(n, j) for j in range(k, n + 1)) / 2**n
    assert p_value <= 0.05, f"unresolved slowing down: n={n}, k={k}, p={p_value}"
    assert np.mean(energy) < 0.99, "wall output lost the slowed-down state"
    print(
        f"Slowing down: n={n}, k={k}, sign-test p={p_value:.6g}; "
        f"mean E/E0={np.mean(energy):.9f}",
        flush=True,
    )
    print("PASS: segment convergence, vacuum limit, slowing down, seeded repeat")
    print(f"Fast wall test: {perf_counter() - begin:.3f} s", flush=True)


if __name__ == "__main__":
    main()
