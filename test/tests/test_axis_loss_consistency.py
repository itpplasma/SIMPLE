"""System-level axis-loss regression for issue #370.

Traces the same deterministic particle set near the axis with the symplectic
Euler integrator (integmode=1) and the RK reference in Boozer coordinates
(integmode=0) on the QA test equilibrium. The old near-axis clamp produced
particles lost only in the symplectic run; the strict acceptance here is
zero symplectic-only losses.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

SIMPLE_IN = """\
&config
trace_time = 1d-2
ntimstep = 200
ntestpart = {npart}
sbeg = 0.05d0
contr_pp = -1e10
multharm = {multharm}
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


def run(simple_x: Path, wout: Path, case: str, integmode: int, relerr: str, startmode: int,
        start_src: Path | None, *, full: bool = False) -> Path:
    d = simple_x.parent / "axis_loss_consistency" / ("full" if full else "fast") / case
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    (d / "wout.nc").symlink_to(wout)
    (d / "simple.in").write_text(
        SIMPLE_IN.format(integmode=integmode, relerr=relerr,
                         startmode=startmode, npart=32 if full else 8,
                         multharm=5 if full else 1))
    if start_src is not None:
        shutil.copy2(start_src, d / "start.dat")
    res = subprocess.run([str(simple_x)], cwd=d, capture_output=True,
                         text=True, timeout=1200)
    if res.returncode != 0:
        print(res.stdout[-2000:])
        print(res.stderr[-2000:])
        sys.exit(f"{case}: simple.x failed")
    return d


def loss_times(d: Path) -> np.ndarray:
    return np.loadtxt(d / "times_lost.dat")[:, 1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("simple_x", type=Path)
    parser.add_argument("wout", type=Path)
    parser.add_argument("--full", action="store_true")
    args = parser.parse_args()
    simple_x = args.simple_x.resolve()
    wout = args.wout.resolve()
    if not simple_x.is_file() or not wout.is_file():
        sys.exit(f"missing test input: {simple_x} or {wout}")

    # Symplectic run generates the deterministic start.dat; RK reuses it.
    d_sympl = run(simple_x, wout, "sympl", integmode=1, relerr="1d-13", startmode=1,
                  start_src=None, full=args.full)
    d_rk = run(simple_x, wout, "rk", integmode=0, relerr="1d-12", startmode=2,
               start_src=d_sympl / "start.dat", full=args.full)

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
