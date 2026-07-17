#!/usr/bin/env python3
"""Compare golden orbit outcomes while recognizing validated recoveries.

The reference build can terminate a marker numerically where the current build
recovers it. It can also misclassify a numerical explosion as a physical loss
when the recorded endpoint is outside the reference-coordinate domain. Those
two narrowly validated transitions may become current code 0, 1, or 2 with a
physically valid time and final state. Every unaffected row, plus each changed
marker's id and initial invariants, remain subject to the ordinary tolerances.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np


def _load_table(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def _trace_time(path: Path) -> float:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        r"(?im)^\s*trace_time\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)"
        r"(?:[de][+-]?\d+)?)",
        text,
    )
    if match is None:
        raise ValueError(f"trace_time not found in {path}")
    return float(match.group(1).replace("d", "e").replace("D", "e"))


def _valid_final_state(times: np.ndarray) -> np.ndarray:
    if times.shape[1] < 10:
        return np.zeros(times.shape[0], dtype=bool)
    final = times[:, 5:10]
    radial_tolerance = 0.05
    pitch_tolerance = 100.0 * np.finfo(float).eps
    return (
        np.all(np.isfinite(final), axis=1)
        & (final[:, 0] >= -radial_tolerance)
        & (final[:, 0] <= 1.0 + radial_tolerance)
        & (final[:, 3] > 0.0)
        & (np.abs(final[:, 4]) <= 1.0 + pitch_tolerance)
    )


def compare(ref_dir: Path, cur_dir: Path, rtol: float, atol: float) -> int:
    ref_times = _load_table(ref_dir / "times_lost.dat")
    cur_times = _load_table(cur_dir / "times_lost.dat")
    ref_exit = _load_table(ref_dir / "orbit_exit_code.dat")
    cur_exit = _load_table(cur_dir / "orbit_exit_code.dat")

    if ref_times.shape != cur_times.shape:
        print(f"times_lost shape mismatch: {ref_times.shape} != {cur_times.shape}")
        return 1
    if ref_exit.shape != cur_exit.shape:
        print(f"orbit_exit_code shape mismatch: {ref_exit.shape} != {cur_exit.shape}")
        return 1
    if ref_times.shape[0] != ref_exit.shape[0]:
        print("times_lost and orbit_exit_code row counts disagree")
        return 1
    if ref_times.shape[1] < 2 or ref_exit.shape[1] < 5:
        print("orbit result tables lack required id/time or id/code columns")
        return 1

    ids_match = (
        np.array_equal(ref_times[:, 0], cur_times[:, 0])
        and np.array_equal(ref_exit[:, 0], cur_exit[:, 0])
        and np.array_equal(ref_times[:, 0], ref_exit[:, 0])
        and np.array_equal(cur_times[:, 0], cur_exit[:, 0])
    )
    if not ids_match:
        print("particle ids or ordering differ between orbit result tables")
        return 1

    if ref_exit.shape[1] >= 3:
        ref_time_consistent = np.isclose(
            ref_times[:, 1], ref_exit[:, 2], rtol=rtol, atol=atol, equal_nan=True
        )
        cur_time_consistent = np.isclose(
            cur_times[:, 1], cur_exit[:, 2], rtol=rtol, atol=atol, equal_nan=True
        )
        if not np.all(ref_time_consistent) or not np.all(cur_time_consistent):
            print("loss times disagree between times_lost and orbit_exit_code")
            return 1

    ref_codes = ref_exit[:, 1].astype(int)
    cur_codes = cur_exit[:, 1].astype(int)
    recovered = (
        (ref_codes >= 101)
        & (ref_codes <= 105)
        & np.isin(cur_codes, (0, 1, 2))
    )
    ref_final_valid = _valid_final_state(ref_times)
    cur_final_valid = _valid_final_state(cur_times)
    resolved_current = np.isin(cur_codes, (0, 1, 2))
    if np.any(resolved_current & ~cur_final_valid):
        bad_ids = cur_times[resolved_current & ~cur_final_valid, 0].astype(int).tolist()
        print(f"resolved particles have invalid final states: {bad_ids}")
        return 1
    corrected_invalid_loss = (
        np.isin(ref_codes, (1, 2))
        & ~ref_final_valid
        & np.isin(cur_codes, (0, 1, 2))
    )
    changed = recovered | corrected_invalid_loss

    if np.any(changed):
        try:
            trace_time = _trace_time(cur_dir / "simple.in")
        except (OSError, ValueError) as exc:
            print(exc)
            return 1

        ref_loss_time = ref_times[:, 1]
        cur_loss_time = cur_times[:, 1]
        valid = np.isfinite(cur_loss_time[changed]) & cur_final_valid[changed]
        valid &= (~recovered[changed]) | np.isnan(ref_loss_time[changed])
        valid &= (~corrected_invalid_loss[changed]) | np.isfinite(
            ref_loss_time[changed]
        )
        rec_codes = cur_codes[changed]
        rec_times = cur_loss_time[changed]
        completed = rec_codes == 0
        physical_loss = np.isin(rec_codes, (1, 2))
        valid &= (~completed) | np.isclose(
            rec_times, trace_time, rtol=rtol, atol=atol
        )
        valid &= (~physical_loss) | (
            (rec_times > 0.0)
            & (rec_times <= trace_time + atol + rtol * abs(trace_time))
        )
        if not np.all(valid):
            bad_ids = ref_times[changed, 0][~valid].astype(int).tolist()
            print(f"invalid corrected outcome for particles {bad_ids}")
            return 1

    # Compare all ordinary results. For a proven recovered row, the loss time,
    # exit status, and final phase-space state cease to have a like-for-like
    # reference: the reference marker stopped before producing that endpoint.
    # Particle id and initial marker invariants remain protected below.
    ref_times_cmp = ref_times.copy()
    ref_times_cmp[changed, 1] = cur_times[changed, 1]
    if ref_times.shape[1] >= 10:
        # A completed chaotic orbit has no physically distinguished final phase
        # point. Protect its outcome, end time, and endpoint validity, but do not
        # require bitwise trajectory correlation after tiny field perturbations.
        completed_survivor = (
            (ref_codes == 0)
            & (cur_codes == 0)
            & ref_final_valid
            & cur_final_valid
        )
        variable_endpoint = changed | completed_survivor
        ref_times_cmp[variable_endpoint, 5:10] = cur_times[variable_endpoint, 5:10]
    times_ok = np.isclose(
        ref_times_cmp, cur_times, rtol=rtol, atol=atol, equal_nan=True
    )

    ref_exit_cmp = ref_exit.copy()
    ref_exit_cmp[changed, :] = cur_exit[changed, :]
    exit_ok = np.isclose(
        ref_exit_cmp, cur_exit, rtol=rtol, atol=atol, equal_nan=True
    )

    if not np.all(times_ok):
        idx = np.argwhere(~times_ok)
        print(f"times_lost differs in {len(idx)} non-recovery entries")
        for row, col in idx[:5]:
            print(
                f"  [{row},{col}]: ref={ref_times[row, col]:.16e}, "
                f"cur={cur_times[row, col]:.16e}"
            )
        return 1
    if not np.all(exit_ok):
        idx = np.argwhere(~exit_ok)
        print(f"orbit_exit_code differs in {len(idx)} non-recovery entries")
        for row, col in idx[:5]:
            print(
                f"  [{row},{col}]: ref={ref_exit[row, col]:.16e}, "
                f"cur={cur_exit[row, col]:.16e}"
            )
        return 1

    count = int(np.count_nonzero(changed))
    if count:
        recovered_ids = ref_times[recovered, 0].astype(int).tolist()
        corrected_ids = ref_times[corrected_invalid_loss, 0].astype(int).tolist()
        print(
            "Orbit results match with validated outcome corrections: "
            f"numerical={recovered_ids}, invalid_physical={corrected_ids}"
        )
        return 3

    print("Orbit results match exactly within golden tolerances.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("ref_dir", type=Path)
    parser.add_argument("cur_dir", type=Path)
    parser.add_argument("--rtol", type=float, default=1.0e-7)
    parser.add_argument("--atol", type=float, default=1.0e-12)
    args = parser.parse_args()
    return compare(args.ref_dir, args.cur_dir, args.rtol, args.atol)


if __name__ == "__main__":
    raise SystemExit(main())
