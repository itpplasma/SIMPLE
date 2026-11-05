#!/usr/bin/env python3
"""Options that a tracing path would silently skip stop at input validation.

Each configuration requests physics that the selected path does not apply:
collisions with SPECTRE fields, a wall with orbit classification, a radial
electric potential with RK recovery, and classification with full orbits.
A volume start whose upper bound lies below the s >= 0.01 clamp would place
every marker outside the requested range. The oracle is the documented contract: simple.x must stop with the matching
message before tracing.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

from simple_driver import link_or_copy

BASE = """
&config
multharm = 3
trace_time = 1d-6
ntestpart = 2
sbeg = 0.5d0
deterministic = .True.
{extra}
/
"""

CASES = (
    ("spectre_collisions",
     "field_input = 'spectre.h5'\nintegmode = 1\nswcoll = .True.",
     "swcoll is not supported for SPECTRE fields"),
    ("wall_classification",
     "netcdffile = 'wout.nc'\nisw_field_type = 2\nwall_input = 'wall.stl'\n"
     "tcut = 5d-7",
     "wall_input is not supported with orbit classification"),
    ("potential_recovery",
     "netcdffile = 'wout.nc'\nisw_field_type = 2\nintegmode = 1\n"
     "radial_electric_potential_slope = 0.1d0",
     "radial_electric_potential_slope needs"),
    ("full_orbit_tcut",
     "netcdffile = 'wout.nc'\nisw_field_type = 2\norbit_model = 7\n"
     "orbit_coord = 1\ntcut = 5d-7",
     "does not support orbit classification"),
    ("volume_start_below_clamp",
     "netcdffile = 'wout.nc'\nisw_field_type = 2\nstartmode = 5\n"
     "num_surf = 1\nsbeg = 0.005d0\ngenerate_start_only = .True.",
     "startmode 5 needs an upper s bound of at least 0.01"),
)


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: test_skipped_physics_rejection.py SIMPLE_X WOUT SPECTRE_H5",
              file=sys.stderr)
        return 2
    simple_x, wout, spectre = (Path(a).resolve() for a in sys.argv[1:])
    work_root = Path(os.environ.get("CTEST_BINARY_DIRECTORY", Path.cwd()))
    failures = []
    for name, extra, message in CASES:
        case_dir = Path(tempfile.mkdtemp(prefix=f"{name}_", dir=work_root))
        (case_dir / "simple.in").write_text(BASE.format(extra=extra))
        link_or_copy(wout, case_dir / "wout.nc")
        link_or_copy(spectre, case_dir / "spectre.h5")
        result = subprocess.run([str(simple_x)], cwd=case_dir, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                timeout=120, check=False)
        if result.returncode == 0 or message not in result.stdout:
            failures.append(f"{name}: expected a stop with '{message}'")
    if failures:
        print("\n".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
