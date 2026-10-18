#!/usr/bin/env python3
"""SPECTRE RK45 single-volume accounting and multi-volume sampling checks."""
import os
import subprocess
import sys
import tempfile

import numpy as np
from spectre_test_utils import axis_volume

NTESTPART = 32
TRACE_TIME = 1.0e-4
SBEG = 0.5
ENERGY_TOL = 1.0e-8        # relative energy drift bound
SAMPLE_TOL = 1.0e-12       # sampler forward-map residual on the interface
PERP_INV_MAX = 1.0e-2      # normalization witness: perp invariant ~ 1/B (Gauss)


def write_input(path, h5, sbeg=SBEG, flux_label=False, generate_only=False):
    lines = [
        "&config",
        f"  trace_time = {TRACE_TIME}",
        f"  sbeg = {sbeg}",
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
    if flux_label:
        lines.insert(-1, "  spectre_sbeg_is_toroidal_flux = .True.")
    if generate_only:
        lines.insert(-1, "  generate_start_only = .True.")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run_simple(binary, workdir, h5, **kwargs):
    write_input(os.path.join(workdir, "simple.in"), h5, **kwargs)
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


def main():
    binary, source = sys.argv[1], sys.argv[2]
    with axis_volume(source) as h5:
        check(binary, h5, source)


def check(binary, h5, source):
    failures = []

    with tempfile.TemporaryDirectory() as work:
        stdout = run_simple(binary, work, h5)
        res = load_times_lost(work)

        lost = res["t_lost"] < TRACE_TIME * (1.0 - 1e-9)
        confined = ~lost
        n_lost = int(lost.sum())
        n_confined = int(confined.sum())

        if n_confined + n_lost != NTESTPART or not np.all(np.isfinite(res["t_lost"])):
            failures.append("single-volume accounting: missing or non-finite marker")
        if n_lost and not np.allclose(res["zend_rho"][lost], 1.0, atol=1e-8):
            failures.append("single-volume loss off outer boundary")
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

        # BDD 3: sampler places markers on the requested surface rho_g = 0.5, and
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

    # Multi-volume surface sampling remains available without orbit tracing.
    samples = []
    for _ in range(2):
        with tempfile.TemporaryDirectory() as work:
            run_simple(binary, work, source, sbeg=2.0, generate_only=True)
            start = np.loadtxt(os.path.join(work, "start.dat"))
            if np.max(np.abs(start[:, 0] - 2.0)) >= SAMPLE_TOL:
                failures.append("multi-volume sampler changed the requested surface")
            samples.append(start)
    if not np.array_equal(*samples):
        failures.append("multi-volume sampling is not deterministic")

    with tempfile.TemporaryDirectory() as work:
        run_simple(binary, work, source, sbeg=0.05, flux_label=True,
                   generate_only=True)
        start = np.loadtxt(os.path.join(work, "start.dat"))
        expected_rho = 0.67072394072143493
        residual = np.max(np.abs(start[:, 0] - expected_rho))
        if not residual < SAMPLE_TOL:
            failures.append(
                f"BDD5: exact flux inverse residual {residual:.3e} >= "
                f"{SAMPLE_TOL:.0e}")
        print(f"BDD5: s=0.05 maps to rho_g={expected_rho:.15f}, "
              f"residual={residual:.3e}")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  " + f)
        sys.exit(1)
    print("\nAll BDD scenarios passed.")


if __name__ == "__main__":
    main()
