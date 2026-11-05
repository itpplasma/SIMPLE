#!/usr/bin/env python3
"""Classification on a logarithmic macrostep grid labels its counts correctly.

The oracle is the run's own loss times: at every row time tau of
confined_fraction.dat, the confined count equals the number of resolved
markers with loss time >= tau (confined markers carry trace_time). Markers
that end unresolved have no loss time and may count either way.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from simple_driver import run_simple_case

CONFIG = """
&config
multharm = 3
contr_pp = -1d10
trace_time = 1d-2
ntimstep = 101
ntestpart = 128
sbeg = 0.6d0
netcdffile = 'wout.nc'
isw_field_type = 2
deterministic = .True.
tcut = 0.5d-2
macrostep_time_grid = 'log'
/
"""


def rows(path: Path) -> list[list[float]]:
    return [[float(x.replace("D", "E")) for x in line.split()]
            for line in path.read_text().splitlines() if line.strip()]


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_classifier_log_macrosteps.py SIMPLE_X WOUT_FILE",
              file=sys.stderr)
        return 2
    case_dir = run_simple_case(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(),
                               "classifier_log_macrosteps", CONFIG, timeout=300)
    loss = [r[1] for r in rows(case_dir / "times_lost.dat")]
    unresolved = sum(1 for t in loss if math.isnan(t))
    failures = []
    for tau, passing, trapped, total in rows(case_dir / "confined_fraction.dat")[1:]:
        confined = round((passing + trapped) * total)
        expected = sum(1 for t in loss
                       if not math.isnan(t) and (t < 0 or t >= tau * (1 - 1e-12)))
        if not expected <= confined <= expected + unresolved:
            failures.append(f"t = {tau:.6e} s: {confined} confined, loss times give {expected}")
    if failures:
        print("\n".join(failures[:5]))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
