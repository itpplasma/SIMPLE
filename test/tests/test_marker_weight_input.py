#!/usr/bin/env python3
"""simple.x writes marker_weight_file rows to results.nc marker_weight.

Wall power is normalised by the source weights, so a run whose markers carry
unequal source power must be able to pass them through simple.x.  Rows that do
not match ntestpart are rejected instead of being padded or truncated.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

SIMPLE_X, WOUT = (Path(arg).resolve() for arg in sys.argv[1:3])
WORK = Path.cwd() / "marker_weight_input"

CONFIG = """&config
trace_time = 1d-6
ntimstep = 2
ntestpart = 4
contr_pp = -1d10
netcdffile = 'wout.nc'
isw_field_type = 2
deterministic = .True.
startmode = 2
integmode = 1
npoiper2 = 32
facE_al = 1.0d0
output_results_netcdf = .True.
{extra}/
"""

START = np.array([
    [0.5, 0.1, 0.2, 1.0, 0.3],
    [0.5, 1.1, 0.4, 1.0, -0.2],
    [0.4, 2.1, 0.6, 0.57735026918962573, 0.5],
    [0.6, 3.1, 0.8, 0.70710678118654757, -0.6],
])


def run(name, extra, weights=None):
    work = WORK / name
    work.mkdir(parents=True, exist_ok=True)
    (work / "results.nc").unlink(missing_ok=True)
    link = work / "wout.nc"
    if not link.exists():
        link.symlink_to(WOUT)
    np.savetxt(work / "start.dat", START, fmt="%.17e")
    if weights is not None:
        np.savetxt(work / "weights.dat", weights, fmt="%.17e")
    (work / "simple.in").write_text(CONFIG.format(extra=extra), encoding="ascii")
    return work, subprocess.run([str(SIMPLE_X)], cwd=work, capture_output=True,
                                text=True, timeout=300)


def main():
    weights = np.array([2.0, 1.0, 0.5, 0.0])
    work, proc = run("supplied", "marker_weight_file = 'weights.dat'\n", weights)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    with Dataset(work / "results.nc") as ds:
        np.testing.assert_array_equal(ds["marker_weight"][:], weights)
        assert ds["marker_weight"].source == "supplied_relative_weights"

    work, proc = run("default", "")
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    with Dataset(work / "results.nc") as ds:
        np.testing.assert_array_equal(ds["marker_weight"][:], 1.0)
        assert ds["marker_weight"].source.startswith("equal_marker_default;")

    for name, rows in (("short", weights[:3]), ("long", np.append(weights, 1.0))):
        work, proc = run(name, "marker_weight_file = 'weights.dat'\n", rows)
        assert proc.returncode != 0, f"{name}: mismatched weight rows accepted"
        assert not (work / "results.nc").exists(), name
    print("marker_weight_file: supplied, default and mismatched rows checked")


if __name__ == "__main__":
    main()
