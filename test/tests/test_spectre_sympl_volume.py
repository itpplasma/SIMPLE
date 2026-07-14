#!/usr/bin/env python3
"""Behavioural test for SPECTRE per-volume symplectic guiding-center tracing (#439).

Runs simple.x on the committed libneo SPECTRE fixture (G3V3L3Fi) with integ_coords
= 6 (SPECTRE) and integmode > 0, markers seeded mid-volume 2 (rho_g = 1.5). The
alpha orbits leave the volume within ~1e-6 s, so a short trace with fine steps
(high npoiper2) keeps every marker confined while resolving ~1e5 symplectic steps.

Scenarios:
  1. Energy: the total energy (p_abs**2 = H) drifts by less than a bound and is
     non-secular (linear-fit slope consistent with zero) over the trace, for the
     symplectic Euler (integmode 1) and midpoint (integmode 3) schemes.
  2. Cross-scheme: the two symplectic schemes integrate the same canonical
     Hamiltonian, so their (rho_g, theta) trajectories agree over the trace.
  3. RK45 cross-check: the direct chart RK45 (integmode 0) and the canonical
     symplectic (integmode 3) both confine every marker to volume 2 and track the
     same orbit initially. They use different guiding-center formulations for
     SPECTRE (RK45 evaluates finite-difference drifts on the stacked chart; the
     symplectic path uses the exact per-volume Meiss Hamiltonian), so they diverge
     at O(rho*) -- a bound, not a tight equality.

Usage:
    test_spectre_sympl_volume.py <simple.x> <spectre_test.h5>
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
from netCDF4 import Dataset

NTESTPART = 24
NTIMSTEP = 200
NPOIPER2 = 6000000        # fine steps -> ~1e5 microsteps within the confined window
TRACE_TIME = 8.0e-8       # shorter than the ~1e-6 s volume-2 confinement time
SBEG = 1.5                # mid-volume 2 (volumes are [0,1], [1,2], [2,3])
RAN_SEED = 42

ENERGY_DRIFT_MAX = 5.0e-6   # bounded energy error over the trace
ENERGY_SECULAR_MAX = 5.0e-6  # non-secular: linear growth over the trace is small
CROSS_S_TOL = 5.0e-3        # symplectic Euler vs midpoint, rho_g
CROSS_TH_TOL = 5.0e-3       # symplectic Euler vs midpoint, theta
RK45_S_TOL = 3.0e-2         # sympl vs chart RK45 over the first stepping macrostep
RK45_TH_TOL = 3.0e-2
RK45_K = 4                  # short-trace window for the RK45 cross-check


def write_input(path, h5, integmode):
    lines = [
        "&config",
        f"  trace_time = {TRACE_TIME}",
        f"  sbeg = {SBEG}",
        f"  ntestpart = {NTESTPART}",
        f"  ntimstep = {NTIMSTEP}",
        f"  npoiper2 = {NPOIPER2}",
        "  relerr = 1d-12",
        f"  field_input = '{h5}'",
        "  integ_coords = 6",
        f"  integmode = {integmode}",
        "  output_orbits_macrostep = .True.",
        "  deterministic = .True.",
        f"  ran_seed = {RAN_SEED}",
        "/",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run_mode(binary, h5, integmode):
    with tempfile.TemporaryDirectory() as work:
        write_input(os.path.join(work, "simple.in"), h5, integmode)
        proc = subprocess.run([binary, "simple.in"], cwd=work,
                              capture_output=True, text=True, timeout=600)
        if proc.returncode != 0:
            raise RuntimeError(
                f"simple.x integmode={integmode} failed:\n"
                f"{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
        ds = Dataset(os.path.join(work, "orbits.nc"))
        # netCDF variables are stored as (timestep, particle).
        s = np.array(ds["s"][:])
        theta = np.array(ds["theta"][:])
        p_abs = np.array(ds["p_abs"][:])
        ds.close()
    return s, theta, p_abs


def confined_mask(p_abs):
    """Particles finite for the whole trace (never boundary-stopped)."""
    return np.all(np.isfinite(p_abs), axis=0)


def energy_stats(p_abs, mask):
    max_drift = 0.0
    max_secular = 0.0
    idx = np.arange(p_abs.shape[0])
    for i in np.nonzero(mask)[0]:
        e = p_abs[:, i] ** 2
        e0 = e[0]
        max_drift = max(max_drift, np.max(np.abs(e - e0)) / e0)
        slope = np.polyfit(idx, e, 1)[0]
        max_secular = max(max_secular, abs(slope) * len(e) / e0)
    return max_drift, max_secular


def traj_agreement(sa, tha, sb, thb, mask, k):
    ds = 0.0
    dth = 0.0
    for i in np.nonzero(mask)[0]:
        a_s, b_s = sa[:k, i], sb[:k, i]
        a_t, b_t = tha[:k, i], thb[:k, i]
        if np.all(np.isfinite(a_s)) and np.all(np.isfinite(b_s)):
            ds = max(ds, np.max(np.abs(a_s - b_s)))
            dth = max(dth, np.max(np.abs(a_t - b_t)))
    return ds, dth


def main():
    binary, h5 = sys.argv[1], sys.argv[2]
    failures = []

    s3, th3, pa3 = run_mode(binary, h5, 3)
    s1, th1, pa1 = run_mode(binary, h5, 1)
    s0, th0, pa0 = run_mode(binary, h5, 0)

    m3 = confined_mask(pa3)
    m1 = confined_mask(pa1)
    m0 = confined_mask(pa0)
    common = m3 & m1 & m0

    if m3.sum() == 0 or m1.sum() == 0:
        failures.append("no marker stayed confined for the whole trace")

    # Scenario 1: bounded, non-secular energy for both symplectic schemes.
    for tag, pa, mask in (("midpoint(3)", pa3, m3), ("euler(1)", pa1, m1)):
        drift, secular = energy_stats(pa, mask)
        if not drift < ENERGY_DRIFT_MAX:
            failures.append(f"energy {tag}: drift {drift:.3e} >= {ENERGY_DRIFT_MAX:.0e}")
        if not secular < ENERGY_SECULAR_MAX:
            failures.append(
                f"energy {tag}: secular {secular:.3e} >= {ENERGY_SECULAR_MAX:.0e}")
        print(f"energy {tag}: confined={int(mask.sum())}/{NTESTPART} "
              f"drift={drift:.3e} secular={secular:.3e}")

    # Scenario 2: the two symplectic schemes track the same canonical orbit.
    ds13, dth13 = traj_agreement(s1, th1, s3, th3, m1 & m3, NTIMSTEP)
    if not ds13 < CROSS_S_TOL:
        failures.append(f"cross-scheme rho_g: {ds13:.3e} >= {CROSS_S_TOL:.0e}")
    if not dth13 < CROSS_TH_TOL:
        failures.append(f"cross-scheme theta: {dth13:.3e} >= {CROSS_TH_TOL:.0e}")
    print(f"cross-scheme 1-vs-3: max|drho_g|={ds13:.3e} max|dtheta|={dth13:.3e}")

    # Scenario 3: sympl and chart RK45 confine every marker and start consistent.
    if int(m0.sum()) != int(m3.sum()):
        failures.append(
            f"RK45/sympl confinement mismatch: {int(m0.sum())} vs {int(m3.sum())}")
    ds30, dth30 = traj_agreement(s3, th3, s0, th0, common, RK45_K)
    if not ds30 < RK45_S_TOL:
        failures.append(f"RK45 rho_g start: {ds30:.3e} >= {RK45_S_TOL:.0e}")
    if not dth30 < RK45_TH_TOL:
        failures.append(f"RK45 theta start: {dth30:.3e} >= {RK45_TH_TOL:.0e}")
    print(f"RK45 cross-check: confined(rk45)={int(m0.sum())} "
          f"start max|drho_g|={ds30:.3e} max|dtheta|={dth30:.3e}")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  " + f)
        sys.exit(1)
    print("\nAll scenarios passed.")


if __name__ == "__main__":
    main()
