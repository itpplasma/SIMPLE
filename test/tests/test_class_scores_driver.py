#!/usr/bin/env python3
"""Driver coverage for the continuous fast-classifier score output.

The independent numerical oracle for the score formulas is the analytic
``test_continuous_classifier_scores`` unit test. This driver verifies that a
real trace writes those scores with their resolution metadata, while retaining
the legacy margins solely for published-classifier compatibility checks.
"""

from __future__ import annotations

import sys
from pathlib import Path

from simple_driver import run_simple_case

NTESTPART = 64
TOL_PERPINV = 15.0

CONFIG = """
&config
netcdffile = 'wout.nc'
multharm = 3
contr_pp = -1d10
trace_time = 4d-2
ntestpart = 64
nper = 1
npoiper = 4
npoiper2 = 64
ntimstep = 16
num_surf = 1
sbeg = 0.3
nturns = 4
notrace_passing = 0
isw_field_type = 2
integmode = 3
deterministic = .True.
fast_class = .True.
tcut = -1d0
/
"""


def read_rows(path: Path) -> list[list[float]]:
    rows = []
    for line in path.read_text().split("\n"):
        if line.strip():
            rows.append([float(value) for value in line.split()])
    return rows


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_class_scores_driver.py SIMPLE_X WOUT_FILE", file=sys.stderr)
        return 2

    case = run_simple_case(
        Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), "class_scores", CONFIG
    )

    scores = read_rows(case / "class_scores.dat")
    parts = read_rows(case / "class_parts.dat")

    if len(scores) != NTESTPART:
        raise AssertionError(f"class_scores.dat has {len(scores)} rows, want {NTESTPART}")
    if len(parts) != NTESTPART:
        raise AssertionError(f"class_parts.dat has {len(parts)} rows, want {NTESTPART}")

    checked = 0
    jpar_scores = 0
    margins = 0
    rotation_scores = 0
    for score, part in zip(scores, parts):
        (index, jpar_rate, rotation_drift, turns, jpar_samples,
         rotation_samples, tip_count, spread, reference, margin, status,
         trap_par) = score
        if int(index) != int(part[0]):
            raise AssertionError("class_scores and class_parts particle order differ")
        if min(jpar_rate, rotation_drift, turns, spread, reference) < 0.0:
            raise AssertionError(f"particle {index}: negative J_parallel magnitude")
        if int(status) not in (0, 1, 2, 3):
            raise AssertionError(f"particle {index}: bad score status {status}")
        if any(value != int(value) for value in (
            jpar_samples, rotation_samples, tip_count, status
        )):
            raise AssertionError(
                f"particle {index}: non-integral resolution metadata"
            )
        if int(jpar_samples) > max(int(tip_count) - 2, 0):
            raise AssertionError(
                f"particle {index}: too many J_parallel samples"
            )
        expected_rotation_samples = (int(tip_count) - 1) // 2
        if expected_rotation_samples < 2:
            expected_rotation_samples = 0
        if int(rotation_samples) != expected_rotation_samples:
            raise AssertionError(f"particle {index}: rotation sample count mismatch")
        if int(jpar_samples) > 0:
            jpar_scores += 1
        if int(rotation_samples) > 0:
            rotation_scores += 1

        # Status 2 never reaches the margin code at all.
        if int(status) not in (1, 3):
            continue
        jpar, ideal = int(part[3]), int(part[4])
        checked += 1

        expected_jpar = 2 if spread > TOL_PERPINV else 1
        if jpar != expected_jpar:
            raise AssertionError(
                f"particle {index}: J_parallel spread {spread} thresholds to "
                f"{expected_jpar} but the class is {jpar}"
            )
        # Only status 1 forms a monotonicity margin; status 3 was decided by
        # the recurrence test, which produces none.
        if int(status) != 1:
            continue
        margins += 1
        expected_ideal = 2 if margin < 0.0 else 1
        if ideal != expected_ideal:
            raise AssertionError(
                f"particle {index}: topology margin {margin} implies "
                f"{expected_ideal} but the class is {ideal}"
            )

    if checked == 0 or jpar_scores == 0 or rotation_scores == 0:
        raise AssertionError(
            "the trace did not resolve both continuous scores and a legacy class"
        )
    # Margin coverage is opportunistic: status 1 needs an orbit the recurrence
    # test calls one-line with a non-empty monotonicity interval, which this
    # field and surface do not reliably produce. The J_parallel comparison is
    # the load-bearing one here.
    print(
        f"continuous J_parallel score on {jpar_scores} orbits; "
        f"rotation drift on {rotation_scores}; legacy threshold agreement on "
        f"{checked}; topology margin exercised on {margins}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
