#!/usr/bin/env python3
"""Driver coverage for standalone fast classification.

fast_class documents itself as "quit immediately after fast classification".
Standalone (tcut <= 0, class_plot = .False.) it is the only configuration that
classifies orbits without the Minkowski fractal cut, which fires at
kt == ntcut, and without disabling the early exit, which requires
.not. class_plot.

The oracle is the classification output itself: a classifying run owes a
class_parts.dat with one row per test particle, and a non-classifying run owes
none. Both expectations are independent of how the dispatch predicate is
spelled.
"""

from __future__ import annotations

import sys
from pathlib import Path

from simple_driver import assert_file_absent, assert_row_count, run_simple_case

NTESTPART = 4

COMMON = """
&config
multharm = 3
contr_pp = -1d10
trace_time = 1d-4
ntestpart = 4
nper = 1
npoiper = 4
npoiper2 = 32
ntimstep = 4
num_surf = 1
sbeg = 0.3
notrace_passing = 1
netcdffile = 'wout.nc'
isw_field_type = 2
integmode = 3
deterministic = .True.
/
"""

# (name, extra namelist lines, classification expected)
CASES = (
    ("fast_class_standalone", "fast_class = .True.\ntcut = -1d0", True),
    ("fast_class_with_cut", "fast_class = .True.\ntcut = 5d-5", True),
    ("fast_class_with_plot", "fast_class = .True.\nclass_plot = .True.", True),
    ("fast_class_absent", "tcut = -1d0", False),
)


def config_with(additions: str) -> str:
    return COMMON.replace(
        "deterministic = .True.",
        f"deterministic = .True.\n{additions}",
    )


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_fast_class_driver.py SIMPLE_X WOUT_FILE", file=sys.stderr)
        return 2

    simple_x = Path(sys.argv[1]).resolve()
    wout = Path(sys.argv[2]).resolve()

    for name, additions, classifies in CASES:
        case_dir = run_simple_case(simple_x, wout, name, config_with(additions))
        if classifies:
            assert_row_count(case_dir, "class_parts.dat", NTESTPART)
        else:
            assert_file_absent(case_dir, "class_parts.dat")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
