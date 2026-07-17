#!/usr/bin/env python3
"""Compare golden orbit outcomes while recognizing recovered markers.

The reference build can terminate a marker numerically where the current build
recovers it.  That one transition is intentional: reference exit code 101--105
and a NaN loss time may become current code 0, 1, or 2 with a physically valid
time and final state. Every unaffected row, plus the recovered marker's id and
initial invariants, remain subject to the ordinary golden tolerances.
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
    if ref_times.shape[1] < 2 or ref_exit.shape[1] < 2:
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

    if np.any(recovered):
        try:
            trace_time = _trace_time(cur_dir / "simple.in")
        except (OSError, ValueError) as exc:
            print(exc)
            return 1

        ref_loss_time = ref_times[:, 1]
        cur_loss_time = cur_times[:, 1]
        valid = np.isnan(ref_loss_time[recovered]) & np.isfinite(
            cur_loss_time[recovered]
        )
        rec_codes = cur_codes[recovered]
        rec_times = cur_loss_time[recovered]
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
            bad_ids = ref_times[recovered, 0][~valid].astype(int).tolist()
            print(f"invalid numerical-recovery outcome for particles {bad_ids}")
            return 1

    # Compare all ordinary results. For a proven recovered row, the loss time,
    # exit status, and final phase-space state cease to have a like-for-like
    # reference: the reference marker stopped before producing that endpoint.
    # Particle id and initial marker invariants remain protected below.
    ref_times_cmp = ref_times.copy()
    ref_times_cmp[recovered, 1] = cur_times[recovered, 1]
    if ref_times.shape[1] >= 10:
        ref_times_cmp[recovered, 5:10] = cur_times[recovered, 5:10]
    times_ok = np.isclose(
        ref_times_cmp, cur_times, rtol=rtol, atol=atol, equal_nan=True
    )

    ref_exit_cmp = ref_exit.copy()
    ref_exit_cmp[recovered, :] = cur_exit[recovered, :]
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

    count = int(np.count_nonzero(recovered))
    if count:
        ids = ref_times[recovered, 0].astype(int).tolist()
        print(f"Orbit results match with {count} numerical recoveries: {ids}")
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
