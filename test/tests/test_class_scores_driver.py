#!/usr/bin/env python3
"""Driver coverage for the continuous classifier margins.

The classifier reduces two continuous quantities to integer codes: the
J_parallel spread is thresholded against tol_perpinv, and the ideal-orbit
class is the sign of the monotonicity margin. class_scores.dat exposes those
quantities so a smooth confinement score can be built without re-tracing.

The oracle is the integer classification SIMPLE already writes: thresholding
the continuous score must reproduce the discrete class it was derived from.
That is a property of the pair of files, independent of either one's values.
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
    margins = 0
    excursions = 0
    for score, part in zip(scores, parts):
        (index, spread, reference, margin, status, radial_spread,
         tip_count, trap_par) = score
        if int(index) != int(part[0]):
            raise AssertionError("class_scores and class_parts particle order differ")
        if spread < 0.0 or reference < 0.0:
            raise AssertionError(f"particle {index}: negative J_parallel magnitude")
        if int(status) not in (0, 1, 2, 3):
            raise AssertionError(f"particle {index}: bad score status {status}")
        # The tip excursion needs no classification, only two tips.
        if radial_spread < 0.0:
            raise AssertionError(f"particle {index}: negative radial excursion")
        if int(tip_count) >= 2 and radial_spread <= 0.0:
            raise AssertionError(
                f"particle {index}: {int(tip_count)} tips but zero excursion"
            )
        if int(tip_count) < 2 and radial_spread != 0.0:
            raise AssertionError(
                f"particle {index}: excursion without two tips"
            )
        if int(tip_count) >= 2:
            excursions += 1

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

    if checked == 0:
        raise AssertionError(
            "no orbit resolved over nturns, so the comparison was vacuous; "
            "raise trace_time or lower nturns in this case"
        )
    # Margin coverage is opportunistic: status 1 needs an orbit the recurrence
    # test calls one-line with a non-empty monotonicity interval, which this
    # field and surface do not reliably produce. The J_parallel comparison is
    # the load-bearing one here.
    if excursions == 0:
        raise AssertionError("no orbit collected two tips, so the excursion is untested")
    print(
        f"J_parallel spread agrees with the class on {checked} orbits; "
        f"topology margin exercised on {margins}; "
        f"radial excursion on {excursions}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
