#!/usr/bin/env python3
"""Reduced ensembles; the unchanged 192-marker test remains the full-tier oracle."""

from concurrent.futures import ProcessPoolExecutor
from math import comb
from pathlib import Path
import sys
from time import perf_counter

from netCDF4 import Dataset
import numpy as np

import test_collisional_stl_wall as full


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


def run(args):
    base, simple, case = args
    name, ensemble, collisions, intervals = case
    start = perf_counter()
    result = full.run_case(
        base,
        name,
        simple,
        ensemble,
        collisions=collisions,
        density=0.0 if name == "dilute" else 2.5e20,
        npoiper2=intervals,
    )
    print(f"{name}: executable {perf_counter() - start:.3f} s", flush=True)
    return result


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
    # Stored collisionless impacts in the 32-marker subset. Running their fine
    # replay independently avoids a second wave of expensive field setup.
    fine_indices = np.array([43, 98, 117, 166, 172])
    cases = (
        ("collisionless", starts[ref_indices], False, 64),
        ("dilute", starts[ref_indices], True, 64),
        ("collisional", starts[coll_indices], True, 64),
        ("repeat", starts[coll_indices], True, 64),
        ("collisionless_fine", starts[fine_indices], False, 128),
    )

    # Each subprocess uses one thread. CTest reserves four CPUs and 5 GiB for
    # the five independent cases (including their overlapping initialization).
    with ProcessPoolExecutor(max_workers=5) as pool:
        reference, dilute, collisional, repeat, fine = list(
            pool.map(run, ((base, simple, case) for case in cases))
        )
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
