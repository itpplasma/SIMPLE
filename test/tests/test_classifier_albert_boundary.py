#!/usr/bin/env python3
"""Classification in Albert coordinates counts a marker leaving s = 1 as lost.

In Albert coordinates s = 1 is not a coordinate surface, so the symplectic
step does not stop there. Plain tracing checks the reference radius after
every step; the classifier did not, so this QA marker, lost at about
6.26e-5 s in plain tracing, ran on outside the plasma until it left the
integration domain: lost about 1e-6 s late, or in longer runs a numerical
failure (exit code 104). The loss times must agree to well within that.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from simple_driver import run_simple_case

START = (
    "5.0000000000000000e-01  3.9872478396860456e+00  7.9996477016781109e-01  "
    "1.0000000000000000e+00  1.5227676977936677e-01\n"
)
CONFIG = """
&config
multharm = 5
ns_s = 5
ns_tp = 5
contr_pp = -1d10
trace_time = {trace_time}
sbeg = 0.5d0
ntestpart = 1
npoiper2 = 256
ntimstep = 100
integmode = 1
relerr = 1d-13
deterministic = .True.
startmode = 2
isw_field_type = 4
netcdffile = 'wout.nc'
tcut = {tcut}
/
"""
ORBIT_EXIT_LCFS = 1


def run(simple_x: Path, wout: Path, name: str, tcut: str):
    config = CONFIG.format(trace_time="2d-4", tcut=tcut)
    case = run_simple_case(simple_x, wout, name, config, timeout=300, start=START)
    code = int((case / "orbit_exit_code.dat").read_text().split()[1])
    lost = float((case / "times_lost.dat").read_text().split()[1].replace("D", "E"))
    return code, lost


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_classifier_albert_boundary.py SIMPLE_X WOUT_FILE", file=sys.stderr)
        return 2
    simple_x, wout = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    failures = []
    plain_code, plain_lost = run(simple_x, wout, "albert_plain", "-1d0")
    class_code, class_lost = run(simple_x, wout, "albert_class", "1d-4")
    if plain_code != ORBIT_EXIT_LCFS:
        failures.append(f"plain tracing: exit code {plain_code}, expected a loss")
    if class_code != ORBIT_EXIT_LCFS:
        failures.append(f"classification: exit code {class_code}, expected a loss")
    elif not math.isclose(class_lost, plain_lost, rel_tol=2e-3):
        failures.append(f"loss times differ: classification {class_lost}, plain {plain_lost}")
    if failures:
        print("\n".join(failures))
        return 1
    print(f"test_classifier_albert_boundary PASSED (lost at {class_lost:.4e} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
