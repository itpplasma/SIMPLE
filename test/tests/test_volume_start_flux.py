#!/usr/bin/env python3
"""Volume starts (startmode = 5) draw s uniformly in the requested range.

startmode 5 with num_surf = 1 documents s uniform in [0.01, sbeg(1)] for every
coordinate set. The oracle is that distribution: all start.dat radii (s in
VMEC reference coordinates) lie in the range and their mean is its midpoint.
"""

from __future__ import annotations

import sys
from pathlib import Path

from simple_driver import run_simple_case

NTEST = 2000
SBEG = 0.5
CONFIG = """
&config
multharm = 3
ntestpart = {ntest}
startmode = 5
num_surf = 1
sbeg = {sbeg}
netcdffile = 'wout.nc'
isw_field_type = {field}
deterministic = .True.
generate_start_only = .True.
/
"""


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_volume_start_flux.py SIMPLE_X WOUT_FILE", file=sys.stderr)
        return 2
    simple_x, wout = (Path(a).resolve() for a in sys.argv[1:])
    low, high = 0.01, SBEG
    tolerance = 4 * (high - low) / (12 * NTEST) ** 0.5
    failures = []
    for name, field in (("boozer", 2), ("meiss", 3), ("albert", 4)):
        case_dir = run_simple_case(
            simple_x, wout, f"volume_start_{name}",
            CONFIG.format(ntest=NTEST, sbeg=SBEG, field=field), timeout=300)
        s = [float(line.split()[0].replace("D", "E"))
             for line in (case_dir / "start.dat").read_text().splitlines() if line.strip()]
        mean = sum(s) / len(s)
        if min(s) < low * (1 - 1e-9) or max(s) > high * (1 + 1e-9) or \
                abs(mean - 0.5 * (low + high)) > tolerance:
            failures.append(f"{name}: s in [{min(s):.4f}, {max(s):.4f}], mean {mean:.4f}")
    if failures:
        print("\n".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
