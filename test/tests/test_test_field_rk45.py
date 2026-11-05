#!/usr/bin/env python3
"""RK45 and a symplectic scheme trace the same TEST-field guiding centres.

The analytic TEST field enters RK45 through magfie_test and the symplectic
schemes through field_can_test. Both describe one guiding-centre model, so
over a short time the end radii of the same markers must agree to the
integrators' accuracy (relative 1e-3 here).
"""

from __future__ import annotations

import sys
from pathlib import Path

from simple_driver import run_simple_case

CONFIG = """
&config
isw_field_type = -1
ntestpart = 16
ntimstep = 201
trace_time = 2d-1
sbeg = 0.3
deterministic = .true.
contr_pp = -1d10
integmode = {integmode}
relerr = 1d-10
npoiper2 = 8192
/
"""


def end_states(case_dir: Path) -> list[list[float]]:
    rows = []
    for line in (case_dir / "times_lost.dat").read_text().splitlines():
        fields = [float(x.replace("D", "E")) for x in line.split()]
        if fields:
            rows.append(fields)
    return rows


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_test_field_rk45.py SIMPLE_X", file=sys.stderr)
        return 2
    simple_x = Path(sys.argv[1]).resolve()
    unused = Path(__file__).resolve()  # the TEST field reads no input file
    rk = end_states(run_simple_case(simple_x, unused, "test_field_rk45",
                                    CONFIG.format(integmode=0)))
    gauss = end_states(run_simple_case(simple_x, unused, "test_field_gauss2",
                                       CONFIG.format(integmode=5)))
    failures = []
    for a, b in zip(rk, gauss):
        if abs(a[5] - b[5]) > 1e-3 * abs(b[5]):
            failures.append(f"marker {int(a[0])}: RK45 r = {a[5]:.6f}, Gauss r = {b[5]:.6f}")
    if failures:
        print("\n".join(failures[:5]))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
