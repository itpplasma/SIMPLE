#!/usr/bin/env python3
"""Classification counts confined markers independently of the thread count.

Passing markers skipped by contr_pp and traced markers update the same
confined-fraction counters from different threads. The oracle is the
single-thread run: confined_fraction.dat must be identical with many threads.
Thousands of short markers make lost updates likely when the counters race.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

from simple_driver import link_or_copy

CONFIG = """
&config
multharm = 3
contr_pp = -1d0
trace_time = 2d-6
ntimstep = 3
ntestpart = 20000
sbeg = 0.5d0
netcdffile = 'wout.nc'
isw_field_type = 2
deterministic = .True.
tcut = 1d-6
/
"""


def run(simple_x: Path, wout: Path, threads: int) -> str:
    work_root = Path(os.environ.get("CTEST_BINARY_DIRECTORY", Path.cwd()))
    case_dir = Path(tempfile.mkdtemp(prefix=f"classifier_threads_{threads}_",
                                     dir=work_root))
    (case_dir / "simple.in").write_text(CONFIG)
    link_or_copy(wout, case_dir / "wout.nc")
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE")
    result = subprocess.run([str(simple_x)], cwd=case_dir, env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=300, check=False)
    if result.returncode != 0:
        raise AssertionError(f"simple.x failed with {threads} threads\n{result.stdout}")
    return (case_dir / "confined_fraction.dat").read_text()


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_classifier_thread_counts.py SIMPLE_X WOUT_FILE",
              file=sys.stderr)
        return 2
    simple_x = Path(sys.argv[1]).resolve()
    wout = Path(sys.argv[2]).resolve()
    serial = run(simple_x, wout, 1)
    threads = max(16, os.cpu_count() or 1)
    parallel = run(simple_x, wout, threads)
    if parallel != serial:
        print(f"confined fractions differ between 1 and {threads} threads")
        print("1 thread:\n" + serial + f"{threads} threads:\n" + parallel)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
