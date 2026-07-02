#!/usr/bin/env python3
"""System-level axis-loss regression for issue #370.

Traces the same deterministic particle set near the axis with the symplectic
Euler integrator (integmode=1) and the RK reference in Boozer coordinates
(integmode=0) on the QA test equilibrium. The old near-axis clamp produced
particles lost only in the symplectic run; the strict acceptance here is
zero symplectic-only losses.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
BUILD_DIR = SCRIPT_DIR.parent.parent / "build"
SIMPLE_X = BUILD_DIR / "simple.x"
WOUT = SCRIPT_DIR.parent / "test_data" / "wout.nc"
WORKDIR = BUILD_DIR / "axis_loss_consistency"

SIMPLE_IN = """\
&config
trace_time = 1d-2
ntimstep = 200
ntestpart = 32
sbeg = 0.05d0
contr_pp = -1e10
netcdffile = 'wout.nc'
isw_field_type = 2
deterministic = .True.
startmode = {startmode}
integmode = {integmode}
relerr = {relerr}
npoiper2 = 256
facE_al = 1.0d0
/
"""


def run(case: str, integmode: int, relerr: str, startmode: int,
        start_src: Path | None) -> Path:
    d = WORKDIR / case
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    (d / "wout.nc").symlink_to(WOUT)
    (d / "simple.in").write_text(
        SIMPLE_IN.format(integmode=integmode, relerr=relerr,
                         startmode=startmode))
    if start_src is not None:
        shutil.copy2(start_src, d / "start.dat")
    res = subprocess.run([str(SIMPLE_X)], cwd=d, capture_output=True,
                         text=True, timeout=1200)
    if res.returncode != 0:
        print(res.stdout[-2000:])
        print(res.stderr[-2000:])
        sys.exit(f"{case}: simple.x failed")
    return d


def loss_times(d: Path) -> np.ndarray:
    return np.loadtxt(d / "times_lost.dat")[:, 1]


def main() -> None:
    if not SIMPLE_X.exists() or not WOUT.exists():
        sys.exit("missing simple.x or wout.nc")

    # Symplectic run generates the deterministic start.dat; RK reuses it.
    d_sympl = run("sympl", integmode=1, relerr="1d-13", startmode=1,
                  start_src=None)
    d_rk = run("rk", integmode=0, relerr="1d-12", startmode=2,
               start_src=d_sympl / "start.dat")

    trace_time = 1e-2
    tl_s = loss_times(d_sympl)
    tl_r = loss_times(d_rk)
    lost_s = (tl_s > 0) & (tl_s < trace_time)
    lost_r = (tl_r > 0) & (tl_r < trace_time)
    sympl_only = lost_s & ~lost_r
    print(f"lost sympl: {lost_s.sum()}, lost RK: {lost_r.sum()}, "
          f"sympl-only: {sympl_only.sum()}, rk-only: {(lost_r & ~lost_s).sum()}")
    if sympl_only.sum() != 0:
        idx = np.nonzero(sympl_only)[0] + 1
        sys.exit(f"FAIL: particles lost only in the symplectic run: {idx}")
    print("test_axis_loss_consistency PASSED")


if __name__ == "__main__":
    main()
