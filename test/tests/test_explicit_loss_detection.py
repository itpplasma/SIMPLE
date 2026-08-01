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

# Two launch surfaces, because the explicit drivers have to re-impose two
# separate chart invariants that the symplectic schemes get from their Newton
# solve.
#
#   0.75  the paper's launch surface, far enough out that a short trace loses
#         markers, which is what exercises the boundary check.
#   0.05  near the axis, which exercises the axis chart switch
#         (r < 0 -> (|r|, theta + pi), issue #370). Integmode 21 records two
#         r_negative events on this configuration, so the path is live; without
#         the switch those markers carry a negative radius into the next
#         substep and the field is evaluated outside the spline domain.
#
#         Be clear about what this surface does and does not test. It runs the
#         chart-switch path and asserts the run completes without faults or
#         escapes, which is a smoke check. It does NOT discriminate whether the
#         switch is applied: removing it changes one marker of 64 by 4.8 in the
#         final state, while two CORRECT methods already differ by up to 5.6 on
#         the same markers, because any trace long enough for an orbit to reach
#         the axis is long enough for the orbits to have decorrelated. No
#         cross-method tolerance separates those two numbers. The switch itself
#         was verified by comparing per-marker final states with it on and off
#         against an otherwise identical build; that evidence lives in the
#         commit, not here. A sensitive regression would need an axis-crossing
#         orbit in the analytic test field, where the trajectory is not chaotic
#         -- worth adding, not attempted here.
#
# contr_pp is forced far negative so no particle is skipped as
# certainly-confined without ever being integrated.
LAUNCH_SURFACES = [0.75, 0.05]
SIMPLE_IN = """\
&config
trace_time = {trace}d0
ntimstep = 100
ntestpart = {npart}
sbeg = {sbeg}d0
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
#
# rtol 1e-6 for the adaptive modes. What is under test is whether a crossing of
# the last closed flux surface is noticed, not how accurately the orbit is
# integrated, and at 1e-6 all five explicit modes already agree with the
# symplectic reference on the loss count to within one marker out of 1000. A
# tighter tolerance buys nothing here and costs a great deal: Gauss-Radau at
# 1e-8 alone took about half an hour of the suite's time. Integmode 24 is
# fixed-step, so its tolerance is inert.
REFERENCE = ("sympl", 3, "1d-13")
EXPLICIT = [
    ("radau", 20, "1d-6"),
    ("gbs", 21, "1d-6"),
    ("cashkarp", 22, "1d-6"),
    ("tdrk", 24, "1d-13"),
    ("tdrk_adaptive", 25, "1d-6"),
]


def run(case, integmode, relerr, startmode, start_src, sbeg):
    d = WORKDIR / case
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    (d / "wout.nc").symlink_to(WOUT)
    (d / "simple.in").write_text(SIMPLE_IN.format(
        trace=TRACE, npart=NPART, integmode=integmode, relerr=relerr,
        startmode=startmode, sbeg=sbeg))
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


def check_surface(sbeg, failures):
    """Run the reference and every explicit mode from one launch surface."""
    tag = f"s{sbeg:g}".replace(".", "")
    name, mode, relerr = REFERENCE
    d_ref = run(f"{tag}_{name}", mode, relerr, startmode=1, start_src=None,
                sbeg=sbeg)
    tl_ref, s_ref = read(d_ref)
    _, lost_ref, _ = classify(tl_ref)
    n_lost_ref = int(np.sum(lost_ref))
    print(f"  sbeg={sbeg:<5g} {name:14s} lost {n_lost_ref:3d}/{NPART}  "
          f"max final s {np.nanmax(s_ref):.4f}")

    # Monte-Carlo tolerance on a loss count of N markers: loose enough that a
    # genuine resolution difference between correct methods passes, far tighter
    # than the 3x gap the missing domain check produced.
    tol = max(4, int(round(3.0 * np.sqrt(max(n_lost_ref, 1)))))

    for name, mode, relerr in EXPLICIT:
        d = run(f"{tag}_{name}", mode, relerr, startmode=2,
                start_src=d_ref / "start.dat", sbeg=sbeg)
        tlost, s = read(d)
        confined, lost, faults = classify(tlost)
        n_lost = int(np.sum(lost))
        n_escaped = int(np.sum(confined & (s > 1.0 + 1e-6)))
        worst = float(np.nanmax(s[confined])) if np.any(confined) else 0.0

        print(f"  sbeg={sbeg:<5g} {name:14s} lost {n_lost:3d}/{NPART}  "
              f"max final s of confined {worst:.4f}  "
              f"escaped {n_escaped}  faulted {faults}")

        if n_escaped > 0:
            failures.append(
                f"sbeg={sbeg:g} {name} (integmode {mode}): {n_escaped} "
                f"marker(s) reported confined while outside the plasma, worst "
                f"at s = {worst:.4f}")
        if abs(n_lost - n_lost_ref) > tol:
            failures.append(
                f"sbeg={sbeg:g} {name} (integmode {mode}): {n_lost} losses "
                f"against the symplectic reference's {n_lost_ref}, outside the "
                f"tolerance of {tol}")
        if faults > 0:
            failures.append(
                f"sbeg={sbeg:g} {name} (integmode {mode}): {faults} orbit(s) "
                f"faulted")
    return n_lost_ref


def main():
    if not SIMPLE_X.exists() or not WOUT.exists():
        sys.exit("missing simple.x or wout.nc")

    failures = []
    losses = [check_surface(sbeg, failures) for sbeg in LAUNCH_SURFACES]

    # The outer surface must actually lose markers, or the boundary check has
    # nothing to detect and a pass would mean nothing. The near-axis surface is
    # there for the chart switch and is not required to lose anything.
    if losses[0] == 0:
        sys.exit("no losses from the outer surface: the test is vacuous")

    if failures:
        for line in failures:
            print(f"FAIL: {line}", file=sys.stderr)
        sys.exit(1)
    print("test_explicit_loss_detection: all checks passed")


if __name__ == "__main__":
    main()
