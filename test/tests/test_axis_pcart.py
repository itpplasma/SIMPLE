#!/usr/bin/env python3
"""Behavioral test for #401: the opt-in pseudo-Cartesian near-axis step.

Trace the same deterministic near-axis particle set with the symplectic Euler1
integrator (integmode=1) on the QA test equilibrium, once with the default flux
chart and once with axis_pcart enabled. The flux chart hits near-axis negative-s
(EVT_R_NEGATIVE) and Newton-maxit (EVT_NEWTON1_MAXIT) faults; axis_pcart must
drive both to zero without losing any particle the default run kept.
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
BUILD_DIR = SCRIPT_DIR.parent.parent / "build"
SIMPLE_X = BUILD_DIR / "simple.x"
WOUT = SCRIPT_DIR.parent / "test_data" / "wout.nc"
WORKDIR = BUILD_DIR / "axis_pcart"
TRACE_TIME = 1e-2

SIMPLE_IN = """\
&config
trace_time = 1d-2
ntestpart = 256
sbeg = 0.05d0
contr_pp = -1e10
multharm = 3
netcdffile = 'wout.nc'
isw_field_type = 2
integmode = 1
npoiper2 = 256
deterministic = .True.
startmode = {startmode}
axis_pcart = {axis_pcart}
axis_pcart_smax = 0.01d0
/
"""


def run(case, axis_pcart, startmode, start_src):
    d = WORKDIR / case
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    (d / "wout.nc").symlink_to(WOUT)
    (d / "simple.in").write_text(
        SIMPLE_IN.format(axis_pcart=axis_pcart, startmode=startmode))
    if start_src is not None:
        shutil.copy2(start_src, d / "start.dat")
    res = subprocess.run([str(SIMPLE_X)], cwd=d, capture_output=True,
                         text=True, timeout=600)
    if res.returncode != 0:
        print(res.stdout[-2000:])
        print(res.stderr[-2000:])
        sys.exit(f"{case}: simple.x failed")
    return d, res.stdout


def faults(stdout):
    """Final cumulative near-axis fault counters from the progress output."""
    def last(name):
        m = re.findall(rf"{name}=(\d+)", stdout)
        return int(m[-1]) if m else 0
    return last("newton1_maxit"), last("r_negative")


def lost(d):
    t = np.loadtxt(d / "times_lost.dat")[:, 1]
    return (t > 0) & (t < TRACE_TIME)


def main():
    if not SIMPLE_X.exists() or not WOUT.exists():
        sys.exit("missing simple.x or wout.nc")

    # Default flux chart generates the deterministic start.dat; axis_pcart reuses it.
    d_off, out_off = run("off", ".False.", startmode=1, start_src=None)
    d_on, out_on = run("on", ".True.", startmode=2, start_src=d_off / "start.dat")

    nm_off, rn_off = faults(out_off)
    nm_on, rn_on = faults(out_on)
    print(f"off: newton1_maxit={nm_off} r_negative={rn_off}")
    print(f"on : newton1_maxit={nm_on} r_negative={rn_on}")

    lost_off, lost_on = lost(d_off), lost(d_on)
    on_only = lost_on & ~lost_off
    print(f"lost off={lost_off.sum()} on={lost_on.sum()} on-only={on_only.sum()}")

    if rn_off + nm_off < 1:
        sys.exit("FAIL: default run hit no near-axis fault; test no longer covers it")
    if rn_on != 0 or nm_on != 0:
        sys.exit(f"FAIL: axis_pcart left faults newton1_maxit={nm_on} r_negative={rn_on}")
    if on_only.sum() != 0:
        sys.exit(f"FAIL: axis_pcart lost particles the default kept: "
                 f"{np.nonzero(on_only)[0] + 1}")

    print("test_axis_pcart PASSED")


if __name__ == "__main__":
    main()
