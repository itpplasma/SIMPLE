#!/usr/bin/env python3
"""Markers that need more than max_rk_recoveries RK recoveries are unresolved.

The analytic TEST tokamak needs a few RK recoveries for 512 markers. With the
default limit no marker is flagged; with max_rk_recoveries = 0 every marker
that needed a recovery ends with exit code 106 (recovery storm, #508) and an
undefined loss time.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from simple_driver import run_simple_case

CONFIG = """
&config
multharm = 5
contr_pp = -1d10
trace_time = 2d-2
sbeg = 0.3d0
ntestpart = 512
npoiper2 = 256
ntimstep = 1000
integmode = 1
relerr = 1d-13
facE_al = 1d0
deterministic = .True.
startmode = 1
isw_field_type = -1
{limit}/
"""
RECOVERY_STORM = 106


def run(simple_x: Path, wout: Path, name: str, limit: str):
    case = run_simple_case(simple_x, wout, name, CONFIG.format(limit=limit), timeout=300)
    codes = [int(l.split()[1]) for l in (case / "orbit_exit_code.dat").read_text().splitlines() if l.strip()]
    lost = [l.split()[1] for l in (case / "times_lost.dat").read_text().splitlines() if l.strip()]
    return codes, lost


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_recovery_storm_guard.py SIMPLE_X WOUT_FILE", file=sys.stderr)
        return 2
    simple_x, wout = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    failures = []

    codes, _ = run(simple_x, wout, "storm_default", "")
    if RECOVERY_STORM in codes:
        failures.append("default limit flagged a recovery storm")

    codes, lost = run(simple_x, wout, "storm_zero", "max_rk_recoveries = 0\n")
    storms = [i for i, c in enumerate(codes) if c == RECOVERY_STORM]
    if not storms:
        failures.append("max_rk_recoveries = 0 flagged no marker although the case needs recoveries")
    for i in storms:
        if not math.isnan(float(lost[i].replace("D", "E"))):
            failures.append(f"storm marker {i + 1} has a finite loss time {lost[i]}")
    if failures:
        print("\n".join(failures))
        return 1
    print(f"test_recovery_storm_guard PASSED ({len(storms)} storm markers)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
