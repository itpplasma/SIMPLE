#!/usr/bin/env python3
"""Classifier output must not depend on the OpenMP thread count.

With num_surf > 1 the classifier looks up per-surface Bmin/Bmax for every
particle. It used to store them in the shared params bmin/bmax while the
thread tracing particle 1 read bmax for the class_plot trapped-passing
boundary (fort.20000), so that file changed from run to run.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

from simple_driver import run_simple_case

CONFIG = """
&config
multharm = 3
contr_pp = -1d10
trace_time = 1d-5
ntestpart = 8
nper = 1
npoiper = 4
npoiper2 = 32
ntimstep = 2
sbeg = 0.2, 0.4
num_surf = 2
netcdffile = 'wout.nc'
isw_field_type = 2
deterministic = .True.
class_plot = .True.
fast_class = .True.
notrace_passing = 1
tcut = -1d0
cut_in_per = 0d0
/
"""

# Sorting preserves duplicate rows while ignoring thread scheduling order.
EXACT = (
    "fort.20000",
    "class_parts.dat",
    "orbit_exit_code.dat",
    "times_lost.dat",
    "confined_fraction.dat",
)
UNORDERED = ("fort.10000", "fort.10011")
THREADS = (1, 4, 4, 4)


def read_lines(case_dir: Path, name: str) -> list[str]:
    path = case_dir / name
    if not path.exists():
        raise AssertionError(f"{case_dir.name}: expected {name}")
    return path.read_text(encoding="utf-8").splitlines()


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_class_plot_threads.py SIMPLE_X WOUT_FILE", file=sys.stderr)
        return 2

    simple_x = Path(sys.argv[1]).resolve()
    wout = Path(sys.argv[2]).resolve()

    cases = []
    for i, nthreads in enumerate(THREADS):
        os.environ["OMP_NUM_THREADS"] = str(nthreads)
        cases.append(
            run_simple_case(simple_x, wout, f"class_plot_t{nthreads}_{i}", CONFIG)
        )

    reference = cases[0]
    failures = []
    for case_dir in cases[1:]:
        for name in EXACT:
            if read_lines(case_dir, name) != read_lines(reference, name):
                failures.append(f"{case_dir.name}: {name} differs from serial run")
        for name in UNORDERED:
            if sorted(read_lines(case_dir, name)) != sorted(
                read_lines(reference, name)
            ):
                failures.append(f"{case_dir.name}: {name} rows differ from serial run")

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
