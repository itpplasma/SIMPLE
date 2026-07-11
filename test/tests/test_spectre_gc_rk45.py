#!/usr/bin/env python3
"""Behavioural test for SPECTRE guiding-center RK45 tracing per volume (#438).

Runs simple.x on the committed libneo SPECTRE fixture with integ_coords='spectre'
(mode id 6) and integmode=0, markers started on control surface rho_g = 2, and
checks the acceptance scenarios:

  1. Every marker either completes the trace inside its volume or boundary-stops
     on an interface to |rho_g - integer| < 1e-10; times_lost.dat and
     confined_fraction.dat account for every marker.
  2. Relative total energy (p^2, conserved by the guiding-center flow) drifts by
     less than the RK45-tolerance bound over the trace.
  3. The surface sampler places markers on the requested interface to < 1e-12 and
     reproduces bit-for-bit with a fixed seed.
  4. The run is VMEC-free: no wout file is configured or opened.

Usage:
    test_spectre_gc_rk45.py <simple.x> <spectre_test.h5>
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

NTESTPART = 32
TRACE_TIME = 1.0e-4
SBEG = 2.0
RHO_TOL = 1.0e-10          # boundary-stop location tolerance
ENERGY_TOL = 1.0e-8        # relative energy drift bound
SAMPLE_TOL = 1.0e-12       # sampler forward-map residual on the interface
PERP_INV_MAX = 1.0e-2      # normalization witness: perp invariant ~ 1/B (Gauss)


def write_input(path, h5):
    lines = [
        "&config",
        f"  trace_time = {TRACE_TIME}",
        f"  sbeg = {SBEG}",
        f"  ntestpart = {NTESTPART}",
        "  ntimstep = 100",
        "  npoiper2 = 256",
        "  relerr = 1d-8",
        f"  field_input = '{h5}'",
        "  integ_coords = 6",
        "  integmode = 0",
        "  deterministic = .True.",
        "  ran_seed = 12345",
        "/",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run_simple(binary, workdir, h5):
    write_input(os.path.join(workdir, "simple.in"), h5)
    proc = subprocess.run([binary, "simple.in"], cwd=workdir,
                          capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        raise RuntimeError(f"simple.x failed:\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout


def load_times_lost(workdir):
    a = np.loadtxt(os.path.join(workdir, "times_lost.dat"))
    return {
        "t_lost": a[:, 1],
        "perp_inv": a[:, 4],
        "zend_rho": a[:, 5],
        "zend_p": a[:, 8],
    }


def count_boundary_events(workdir):
    path = os.path.join(workdir, "spectre_boundary_events.dat")
    n = 0
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def main():
    binary, h5 = sys.argv[1], sys.argv[2]
    failures = []

    with tempfile.TemporaryDirectory() as work:
        stdout = run_simple(binary, work, h5)
        res = load_times_lost(work)
        n_events = count_boundary_events(work)

        boundary = res["t_lost"] < TRACE_TIME * (1.0 - 1e-9)
        confined = ~boundary
        n_boundary = int(boundary.sum())
        n_confined = int(confined.sum())

        # BDD 1: accounting -- every marker confined or boundary-stopped.
        if n_confined + n_boundary != NTESTPART:
            failures.append(
                f"BDD1: {n_confined}+{n_boundary} != {NTESTPART} markers")
        if n_boundary != n_events:
            failures.append(
                f"BDD1: {n_boundary} boundary times but {n_events} logged events")
        if n_boundary:
            worst = np.max(np.abs(res["zend_rho"][boundary]
                                  - np.round(res["zend_rho"][boundary])))
            if not worst < RHO_TOL:
                failures.append(
                    f"BDD1: boundary-stop |rho_g-int| = {worst:.3e} >= {RHO_TOL:.0e}")
        else:
            worst = 0.0
        if n_confined:
            tmin = res["t_lost"][confined].min()
            if not tmin >= TRACE_TIME * (1.0 - 1e-6):
                failures.append(
                    f"BDD1: confined marker loss time {tmin:.3e} < trace_time")
        print(f"BDD1: confined={n_confined} boundary={n_boundary} "
              f"max|rho_g-int|={worst:.3e}")

        # confined_fraction.dat closes the account at the final time.
        cf = np.loadtxt(os.path.join(work, "confined_fraction.dat"))
        frac_final = cf[-1, 1] + cf[-1, 2]
        expected = n_confined / NTESTPART
        if abs(frac_final - expected) > 1.0 / NTESTPART * 0.5 + 1e-9:
            failures.append(
                f"BDD1: confined_fraction {frac_final:.4f} != {expected:.4f}")
        if int(cf[-1, 3]) != NTESTPART:
            failures.append(f"BDD1: confined_fraction n = {int(cf[-1, 3])}")

        # BDD 2: energy. p is a constant of the guiding-center flow, so the
        # relative drift of the total energy (~ p^2) stays at RK45 roundoff.
        drift = np.max(np.abs(res["zend_p"] ** 2 - 1.0))
        if not drift < ENERGY_TOL:
            failures.append(f"BDD2: energy drift {drift:.3e} >= {ENERGY_TOL:.0e}")
        print(f"BDD2: max relative energy drift = {drift:.3e}")

        # Normalization witness: perp invariant ~ v_perp^2/B with B in Gauss.
        # Dropping the Tesla->Gauss factor inflates this ~1e4x (red proof).
        perp_max = res["perp_inv"].max()
        if not perp_max < PERP_INV_MAX:
            failures.append(
                f"Norm: perp_inv max {perp_max:.3e} >= {PERP_INV_MAX:.0e} "
                f"(Gaussian B mis-scaled?)")
        print(f"Norm: perp_inv max = {perp_max:.3e} (< {PERP_INV_MAX:.0e})")

        # BDD 3: sampler places markers exactly on interface rho_g = 2, and
        # reproduces bit-for-bit with the fixed seed.
        start = np.loadtxt(os.path.join(work, "start.dat"))
        res_map = np.max(np.abs(start[:, 0] - SBEG))
        if not res_map < SAMPLE_TOL:
            failures.append(
                f"BDD3: sampler |rho_g-{SBEG}| = {res_map:.3e} >= {SAMPLE_TOL:.0e}")
        with tempfile.TemporaryDirectory() as work2:
            run_simple(binary, work2, h5)
            start2 = np.loadtxt(os.path.join(work2, "start.dat"))
            if not np.array_equal(start, start2):
                failures.append("BDD3: sampler not reproducible with fixed seed")
        print(f"BDD3: sampler |rho_g-{SBEG}| = {res_map:.3e}, reproducible")

        # BDD 4: VMEC-free. No wout is configured, none is written, and the log
        # never enters the VMEC initialization phase.
        stray = [f for f in os.listdir(work)
                 if f.startswith("wout") and f.endswith(".nc")]
        if stray:
            failures.append(f"BDD4: stray VMEC file(s) {stray}")
        if "VMEC initialization completed" in stdout:
            failures.append("BDD4: run entered the VMEC initialization phase")
        if "SPECTRE field loading completed" not in stdout:
            failures.append("BDD4: SPECTRE field loading phase missing from log")
        print("BDD4: VMEC-free initialization confirmed")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  " + f)
        sys.exit(1)
    print("\nAll BDD scenarios passed.")


if __name__ == "__main__":
    main()
