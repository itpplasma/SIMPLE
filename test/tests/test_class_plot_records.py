#!/usr/bin/env python3
"""class_plot writes one record per marker and class.

With class_plot, every classified marker appears once in the Minkowski,
J_parallel and topology files (fort.100xx, fort.400xx, fort.500xx). The
oracle is uniqueness: no file may repeat a marker record, and the Minkowski
files together hold at most ntestpart records.
"""

from __future__ import annotations

import sys
from pathlib import Path

from simple_driver import run_simple_case

NTESTPART = 32
CONFIG = f"""
&config
notrace_passing = 1
ntestpart = {NTESTPART}
num_surf = 0
sbeg = 0.2
npoiper2 = 128
netcdffile = 'wout.nc'
isw_field_type = 2
trace_time = 1.5d-2
tcut = 1d-2
class_plot = .True.
cut_in_per = 0d0
fast_class = .True.
deterministic = .True.
/
"""


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_class_plot_records.py SIMPLE_X WOUT_FILE", file=sys.stderr)
        return 2
    case_dir = run_simple_case(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(),
                               "class_plot_records", CONFIG, timeout=300)
    failures = []
    minkowski = 0
    for path in sorted(case_dir.glob("fort.*")):
        unit = int(path.suffix[1:])
        if not (10011 <= unit <= 10022 or 40012 <= unit <= 40032
                or 50012 <= unit <= 50032):
            continue
        rows = [line for line in path.read_text().splitlines() if line.strip()]
        if len(rows) != len(set(rows)):
            failures.append(f"{path.name}: {len(rows)} records, {len(set(rows))} distinct")
        if unit < 20000:
            minkowski += len(rows)
    if minkowski > NTESTPART:
        failures.append(f"{minkowski} Minkowski records for {NTESTPART} markers")
    if failures:
        print("\n".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
