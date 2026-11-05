#!/usr/bin/env python3
"""ANTS starts (startmode = 3) take one marker per input line.

The oracle is the input itself: marker i starts on the flux surface s given in
line i of start_ants.dat, and times_lost.dat reports that start surface in its
fourth column.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

from simple_driver import link_or_copy

SURFACES = (0.25, 0.36, 0.49)

CONFIG = """
&config
multharm = 3
contr_pp = -1d10
trace_time = 1d-6
ntimstep = 2
ntestpart = 3
startmode = 3
special_ants_file = .true.
netcdffile = 'wout.nc'
isw_field_type = 2
deterministic = .True.
/
"""


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_ants_start_rows.py SIMPLE_X WOUT_FILE", file=sys.stderr)
        return 2
    simple_x = Path(sys.argv[1]).resolve()
    wout = Path(sys.argv[2]).resolve()
    work_root = Path(os.environ.get("CTEST_BINARY_DIRECTORY", Path.cwd()))
    case_dir = Path(tempfile.mkdtemp(prefix="ants_start_rows_", dir=work_root))
    (case_dir / "simple.in").write_text(CONFIG)
    (case_dir / "start_ants.dat").write_text("".join(
        f"v_par = {0.2 * (k + 1)} v_perp = 1.0 u = {0.1 * (k + 1)} "
        f"v = {0.3 * (k + 1)} s = {s}\n" for k, s in enumerate(SURFACES)))
    link_or_copy(wout, case_dir / "wout.nc")
    result = subprocess.run([str(simple_x)], cwd=case_dir, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=120, check=False)
    if result.returncode != 0:
        print(result.stdout)
        return 1
    starts = [float(line.split()[3].replace("D", "E"))
              for line in (case_dir / "times_lost.dat").read_text().splitlines()
              if line.strip()]
    if any(abs(got - want) > 1e-12 for got, want in zip(starts, SURFACES)) or \
            len(starts) != len(SURFACES):
        print(f"start surfaces {starts}, expected {list(SURFACES)}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
