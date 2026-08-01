#!/usr/bin/env python3
"""The explicit integmodes must notice when a marker leaves the plasma.

The symplectic schemes test the radius of every Newton iterate, so a crossing
of the last closed flux surface is caught inside the step. The explicit modes
(20 Gauss-Radau, 21 Bulirsch-Stoer, 22 Cash-Karp, 24/25 TDRK) hand a whole
substep to an adaptive integrator in one call and see nothing between its
endpoints. Until they reported a domain exit, a marker that crossed s = 1 was
integrated onwards through extrapolated splines and recorded as confined: at
paper scale that put 44 of 1000 markers outside the plasma, one as far out as
s = 2.8, and the explicit modes reported 76 losses where the symplectic
schemes and RK4/5 both reported about 207.

The oracle is a physical invariant, not a comparison against a recorded run:

    a marker reported as confined is inside the last closed flux surface.

s > 1 is outside the plasma by definition, so a confined marker sitting there
is a contradiction no tolerance can excuse. This needs no reference
trajectory, and one escaped marker is enough to fail it -- which matters,
because a short trace produces few losses and a test that needed many would be
insensitive.

A cross-method check on the loss count is included as a second, weaker signal:
it depends on the symplectic run being right, but it catches systematic
under-detection that the invariant alone would miss if markers happened to
drift back inside before the trace ended.
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
WORKDIR = BUILD_DIR / "explicit_loss_detection"

NPART = 64
TRACE = 1e-2

# sbeg = 0.75 is the paper's launch surface and far enough out that a short
# trace still loses markers; near the axis almost nothing is lost and the test
# would have nothing to detect. contr_pp is forced far negative so no particle
# is skipped as certainly-confined without ever being integrated.
SIMPLE_IN = """\
&config
trace_time = {trace}d0
ntimstep = 100
ntestpart = {npart}
sbeg = 0.75d0
contr_pp = -1d10
netcdffile = 'wout.nc'
isw_field_type = 2
deterministic = .True.
startmode = {startmode}
integmode = {integmode}
relerr = {relerr}
npoiper2 = 32
facE_al = 1.0d0
/
"""

# Reference plus every explicit mode. The reference generates the start file;
# the rest reuse it, so all of them trace the same markers.
REFERENCE = ("sympl", 3, "1d-13")
EXPLICIT = [
    ("radau", 20, "1d-8"),
    ("gbs", 21, "1d-8"),
    ("cashkarp", 22, "1d-8"),
    ("tdrk", 24, "1d-13"),
    ("tdrk_adaptive", 25, "1d-8"),
]


def run(case, integmode, relerr, startmode, start_src):
    d = WORKDIR / case
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    (d / "wout.nc").symlink_to(WOUT)
    (d / "simple.in").write_text(SIMPLE_IN.format(
        trace=TRACE, npart=NPART, integmode=integmode, relerr=relerr,
        startmode=startmode))
    if start_src is not None:
        shutil.copy2(start_src, d / "start.dat")
    res = subprocess.run([str(SIMPLE_X)], cwd=d, capture_output=True,
                         text=True, timeout=1800)
    if res.returncode != 0:
        print(res.stdout[-2000:])
        print(res.stderr[-2000:])
        sys.exit(f"{case}: simple.x failed")
    return d


def read(d):
    """(loss time, final normalised flux) per marker."""
    data = np.loadtxt(d / "times_lost.dat")
    return data[:, 1], data[:, 5]


def classify(tlost):
    """Confined markers reached the end of the trace; lost ones did not.

    SIMPLE writes times_lost = trace_time for a marker that survives rather
    than a sentinel, so the comparison has to be against the trace time. NaN
    marks an integrator fault, which is neither.
    """
    finite = ~np.isnan(tlost)
    confined = finite & (tlost >= 0.999999 * TRACE)
    lost = finite & (tlost > 0) & (tlost < 0.999999 * TRACE)
    return confined, lost, int(np.sum(~finite))


def main():
    if not SIMPLE_X.exists() or not WOUT.exists():
        sys.exit("missing simple.x or wout.nc")

    name, mode, relerr = REFERENCE
    d_ref = run(name, mode, relerr, startmode=1, start_src=None)
    tl_ref, s_ref = read(d_ref)
    conf_ref, lost_ref, faults_ref = classify(tl_ref)
    n_lost_ref = int(np.sum(lost_ref))
    print(f"  {name:14s} lost {n_lost_ref:3d}/{NPART}  "
          f"max final s {np.nanmax(s_ref):.4f}")

    if n_lost_ref == 0:
        sys.exit("reference lost no markers: the test cannot detect anything")

    # Monte-Carlo tolerance on a loss count of N markers, generous enough that
    # a genuine resolution difference between correct methods does not fail it
    # but far tighter than the 3x gap the missing domain check produced.
    tol = max(4, int(round(3.0 * np.sqrt(max(n_lost_ref, 1)))))

    failures = []
    for name, mode, relerr in EXPLICIT:
        d = run(name, mode, relerr, startmode=2,
                start_src=d_ref / "start.dat")
        tlost, s = read(d)
        confined, lost, faults = classify(tlost)
        n_lost = int(np.sum(lost))
        escaped = confined & (s > 1.0 + 1e-6)
        n_escaped = int(np.sum(escaped))
        worst = float(np.nanmax(s[confined])) if np.any(confined) else 0.0

        print(f"  {name:14s} lost {n_lost:3d}/{NPART}  "
              f"max final s of confined markers {worst:.4f}  "
              f"escaped {n_escaped}  faulted {faults}")

        if n_escaped > 0:
            failures.append(
                f"{name} (integmode {mode}): {n_escaped} marker(s) reported "
                f"confined while outside the plasma, worst at s = {worst:.4f}")
        if abs(n_lost - n_lost_ref) > tol:
            failures.append(
                f"{name} (integmode {mode}): {n_lost} losses against the "
                f"symplectic reference's {n_lost_ref}, outside the tolerance "
                f"of {tol}")

    if failures:
        for line in failures:
            print(f"FAIL: {line}", file=sys.stderr)
        sys.exit(1)
    print("test_explicit_loss_detection: all checks passed")


if __name__ == "__main__":
    main()
