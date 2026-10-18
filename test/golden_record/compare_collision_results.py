#!/usr/bin/env python3
"""Integrity contract for the 32-marker collision golden fixture.

Historical correlated kicks are not an independent physical distribution oracle.
This checks current-run accounting and detects repeated survivor momentum kicks;
analytic collision tests and the collisional wall regression supply physics evidence.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import numpy as np

from compare_orbit_results import _load_table, _trace_time, _valid_final_state


def _fixture(path: Path) -> float:
    text = path.read_text(encoding="utf-8")
    # wall_input defaults to ''; wall_units alone does not enable a wall.
    literal = r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\""
    uncommented = re.sub(
        rf"({literal})|![^\n]*", lambda match: match.group(1) or "", text
    )
    for assignment in re.finditer(
        r"(?i)\bwall_input\s*(?:\([^)]*\)\s*)?=\s*", uncommented
    ):
        value = re.match(literal, uncommented[assignment.end() :])
        if value is None or value.group()[1:-1].strip():
            raise ValueError(f"{path}: collision fixture must not configure a wall")
    expected = {"ntestpart": "32", "swcoll": ".true.", "deterministic": ".true."}
    for key, value in expected.items():
        values = re.findall(rf"(?im)^\s*{key}\s*=\s*([^!\n/]+)", text)
        if len(values) != 1 or values[0].strip().lower() != value:
            raise ValueError(f"{path}: collision fixture requires {key}={value}")
    duration = _trace_time(path)
    if not np.isfinite(duration) or duration <= 0:
        raise ValueError("collision fixture requires a positive finite horizon")
    return duration


def _historical_summary(times, exits):
    codes = exits[:, 1]
    energy = times[:, 8] ** 2
    valid = np.isfinite(energy) & np.isin(codes, [0, 1])
    mean = float(energy[valid].mean()) if np.any(valid) else float("nan")
    return (
        int((codes == 1).sum()),
        int((~np.isin(codes, [0, 1])).sum()),
        mean,
    )


def compare(ref_dir: Path, cur_dir: Path, rtol: float, atol: float) -> int:
    try:
        duration = _fixture(cur_dir / "simple.in")
        if not np.isclose(
            _fixture(ref_dir / "simple.in"), duration, rtol=rtol, atol=atol
        ):
            raise ValueError("reference/current collision horizons differ")
        tables = []
        for path in (ref_dir, cur_dir):
            starts = _load_table(path / "start.dat")
            times = _load_table(path / "times_lost.dat")
            exits = _load_table(path / "orbit_exit_code.dat")
            if (
                starts.shape != (32, 5)
                or times.shape != (32, 10)
                or exits.shape != (32, 5)
            ):
                raise ValueError("collision fixture requires all 32 start/outcome rows")
            ids = np.arange(1, 33)
            if not np.array_equal(times[:, 0], ids) or not np.array_equal(
                exits[:, 0], ids
            ):
                raise ValueError("missing, duplicated or reordered marker IDs")
            if not np.isfinite(starts).all() or not np.isfinite(times[:, 2:5]).all():
                raise ValueError("nonfinite starts or initial invariants")
            if np.any(starts[:, 3] != 1) or np.any(np.abs(starts[:, 4]) > 1):
                raise ValueError(
                    "collision fixture requires unit launch momentum and valid pitch"
                )
            if not np.allclose(
                times[:, 1], exits[:, 2], rtol=rtol, atol=atol, equal_nan=True
            ):
                raise ValueError("loss-time and exit-code files disagree")
            tables.append((starts, times, exits))
        ref_starts, ref_times, ref_exit = tables[0]
        starts, times, exits = tables[1]
        if not np.allclose(starts, ref_starts, rtol=rtol, atol=atol):
            raise ValueError("initial markers changed")
        if not np.allclose(times[:, 2:5], ref_times[:, 2:5], rtol=rtol, atol=atol):
            raise ValueError("initial marker invariants changed")
        if not np.allclose(times[:, 3], starts[:, 0], rtol=rtol, atol=atol):
            raise ValueError("reported initial radius disagrees with start.dat")
        if not np.isfinite(times).all() or not np.isfinite(exits).all():
            raise ValueError("nonfinite current outcome")
        codes = exits[:, 1]
        if not np.isin(codes, [0, 1]).all():
            raise ValueError(
                "impossible wall or unaccounted numerical/skipped/core outcome"
            )
        if not _valid_final_state(times).all():
            raise ValueError("invalid current endpoint, momentum or pitch")
        if np.any(np.abs(times[codes == 1, 5] - 1.0) > 0.05):
            raise ValueError("LCFS exit endpoint is not near the boundary")
        if np.any(times[:, 1] <= 0) or np.any(
            times[:, 1] > duration + atol + rtol * duration
        ):
            raise ValueError("current outcome time lies outside trace horizon")
        complete = codes == 0
        if not np.allclose(times[complete, 1], duration, rtol=rtol, atol=atol):
            raise ValueError("completed marker did not reach trace horizon")
        momenta = times[complete, 8]
        # Flat nonzero collision coefficients and equal start momenta in this
        # fixture made repeated per-marker seeds produce duplicated p histories.
        # Use roundoff-scale equality, not the much wider golden trajectory rtol.
        if len(momenta) >= 2:
            ordered = np.sort(momenta)
            threshold = 64 * np.finfo(float).eps * np.maximum(1.0, np.abs(ordered[1:]))
            if np.any(np.diff(ordered) <= threshold):
                raise ValueError(
                    "duplicated survivor momenta: possible repeated collision stream"
                )
    except (OSError, ValueError, IndexError) as exc:
        print(f"Collision integrity contract failed: {exc}")
        return 1
    ref_loss, ref_unresolved, ref_mean = _historical_summary(ref_times, ref_exit)
    cur_loss, cur_unresolved, cur_mean = _historical_summary(times, exits)
    print("Collision integrity contract passed; this is not distribution validation.")
    print(
        f"Descriptive historical movement: losses {ref_loss}/32 -> {cur_loss}/32; "
        f"unresolved {ref_unresolved} -> {cur_unresolved}; mean resolved p^2 "
        f"{ref_mean:.12g} -> {cur_mean:.12g} (delta={cur_mean - ref_mean:+.12g})."
    )
    print(
        "Required physics evidence: test_collision_substeps, test_particle_rng, "
        "and test_collisional_stl_wall; see docs/collision-golden-contract.md."
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ref_dir", type=Path)
    parser.add_argument("cur_dir", type=Path)
    parser.add_argument("--rtol", type=float, default=1.0e-7)
    parser.add_argument("--atol", type=float, default=1.0e-12)
    args = parser.parse_args()
    return compare(args.ref_dir, args.cur_dir, args.rtol, args.atol)


if __name__ == "__main__":
    raise SystemExit(main())
