#!/usr/bin/env python3
"""Classification starts markers where plain tracing starts them.

trace_orbit maps start.dat from reference (VMEC) coordinates into the
integration coordinates; the classifier must do the same for every
canonical coordinate set. The oracle is the plain trace of the same
markers in the same coordinates: a marker the plain trace loses before
the classification cut tcut must be lost (or end unresolved) in the
classifying run, and the classifying run must lose no marker that the
plain trace keeps.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from simple_driver import run_simple_case

TRACE_TIME = 5e-4
TCUT = 2.5e-4

COMMON = f"""
&config
multharm = 3
contr_pp = -1d10
trace_time = {TRACE_TIME}
ntestpart = 32
canonical_grid_nr = 24
canonical_grid_ntheta = 24
canonical_grid_nphi = 24
sbeg = 0.5d0
netcdffile = 'wout.nc'
deterministic = .True.
"""


def loss_times(case_dir: Path) -> list[float]:
    times = []
    for line in (case_dir / "times_lost.dat").read_text().splitlines():
        fields = line.split()
        if len(fields) > 1:
            times.append(float(fields[1].replace("D", "E")))
    return times


def lost(times: list[float], before: float) -> set[int]:
    return {i for i, t in enumerate(times) if 0.0 < t < before}


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_classifier_start_coordinates.py SIMPLE_X WOUT_FILE",
              file=sys.stderr)
        return 2
    simple_x = Path(sys.argv[1]).resolve()
    wout = Path(sys.argv[2]).resolve()

    failures = []
    for name, field_type in (("meiss", 3), ("albert", 4)):
        config = COMMON + f"isw_field_type = {field_type}\n"
        plain = loss_times(run_simple_case(
            simple_x, wout, f"plain_{name}", config + "tcut = -1d0\n/\n",
            timeout=300))
        classified = loss_times(run_simple_case(
            simple_x, wout, f"classified_{name}", config + f"tcut = {TCUT}\n/\n",
            timeout=300))
        unresolved = {i for i, t in enumerate(classified) if math.isnan(t)}
        prompt = lost(plain, TCUT)
        cut_lost = lost(classified, TRACE_TIME * (1.0 - 1e-9))
        if not prompt:
            failures.append(f"{name}: no prompt losses, the oracle is empty")
        if not prompt <= cut_lost | unresolved:
            failures.append(f"{name}: plain-trace losses {sorted(prompt - cut_lost)} "
                            "survive classification")
        if not cut_lost <= lost(plain, TRACE_TIME * (1.0 - 1e-9)):
            failures.append(f"{name}: classification loses markers "
                            f"{sorted(cut_lost - lost(plain, TRACE_TIME))} "
                            "that the plain trace keeps")
    if failures:
        print("\n".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
