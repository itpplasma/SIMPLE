#!/usr/bin/env python3
"""Behavioural test for the SPECTRE symplectic field-line Poincare app.

Runs spectre_poincare.x on the committed libneo G3V3L3Fi fixture and compares
its per-volume sections and rotational-transform profiles against SPECTRE's own
tracer reference (test/tests/data/poincare.npz).

Usage:
    test_spectre_poincare.py <spectre_poincare.x> <spectre_test.h5> <poincare.npz>

The app seeds field lines at theta = 0.2 with s spread uniformly across each
volume, while the reference used its own seeds (theta = 0, and s in [-0.5, 0.85]
for the axis volume). We therefore compare cloud geometry, not point pairs, and
key the precise geometric check on the outermost seed s = 0.85, the only seed
radius shared by both seed sets in every volume.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

# The reference sections floor at ~3e-3 m relative to our tracer (SPECTRE uses a
# different integrator and radial resolution), so 5e-3 m leaves real margin.
NN_MEDIAN_TOL = 5.0e-3
IOTA_EXPAND = 0.10
ENV_OVERSHOOT = 0.02
ENV_REACH = 0.05


def write_namelist(path, h5, nturns, nsteps, nseeds, volumes=None,
                   seed_file=None):
    lines = ["&poincare", f"  spectre_file = '{h5}'", "  zeta0 = 0.0",
             f"  nturns = {nturns}", f"  nsteps_per_period = {nsteps}",
             f"  nseeds_per_volume = {nseeds}"]
    if volumes is not None:
        lines.append("  volumes = " + ", ".join(str(v) for v in volumes))
    if seed_file is not None:
        lines.append(f"  seed_file = '{seed_file}'")
    lines.append("/")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run_app(binary, workdir, h5, **kw):
    nml = os.path.join(workdir, "spectre_poincare.in")
    write_namelist(nml, h5, **kw)
    proc = subprocess.run([binary, nml], cwd=workdir, capture_output=True,
                          text=True, timeout=600)
    if proc.returncode != 0:
        raise RuntimeError(f"app failed:\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout


def load_volume(workdir, lvol):
    a = np.loadtxt(os.path.join(workdir, f"poincare_vol{lvol}.dat"))
    return a  # columns: seed, section, s, theta, R, Z


def load_iota(workdir, lvol):
    vals = []
    with open(os.path.join(workdir, f"iota_vol{lvol}.dat")) as f:
        for line in f:
            if line.startswith("#"):
                continue
            tok = line.split()
            try:
                vals.append(float(tok[1]))
            except (IndexError, ValueError):
                pass  # "terminated" seeds carry no iota
    return np.array(vals)


def median_nearest(app_pts, ref_pts):
    dmin = np.array([np.hypot(ref_pts[:, 0] - p[0],
                              ref_pts[:, 1] - p[1]).min() for p in app_pts])
    return np.median(dmin)


def final_section_point(a, seed):
    m = a[a[:, 0] == seed]
    k = m[:, 1].max()
    row = m[m[:, 1] == k][0]
    return row[2], row[3]  # s, theta


def main():
    binary, h5, npz_path = sys.argv[1], sys.argv[2], sys.argv[3]
    ref = np.load(npz_path)
    nvol = int(ref["mvol"])

    # Single geometric centre for radial ordering across all volumes.
    all_r = np.concatenate([ref[f"sec_R_vol{i}"].ravel() for i in range(nvol)])
    all_z = np.concatenate([ref[f"sec_Z_vol{i}"].ravel() for i in range(nvol)])
    r0, z0 = all_r.mean(), all_z.mean()

    failures = []
    with tempfile.TemporaryDirectory() as work:
        stdout = run_app(binary, work, h5, nturns=60, nsteps=256, nseeds=4)

        # BDD 1: confinement -- no seed leaves its volume (B.n = 0 witness).
        if "0 seed(s) terminated" not in stdout:
            failures.append(f"BDD1: a seed terminated:\n{stdout}")
        smax = 0.0
        for lvol in range(1, nvol + 1):
            smax = max(smax, np.abs(load_volume(work, lvol)[:, 2]).max())
        if not smax < 1.0:
            failures.append(f"BDD1: max|s| = {smax:.6f} reached the interface")
        print(f"BDD1: max|s| = {smax:.5f}  (no terminations)")

        for lvol in range(1, nvol + 1):
            i = lvol - 1
            a = load_volume(work, lvol)
            ref_pts = np.column_stack([ref[f"sec_R_vol{i}"].ravel(),
                                       ref[f"sec_Z_vol{i}"].ravel()])
            ref_rad = np.hypot(ref_pts[:, 0] - r0, ref_pts[:, 1] - z0)
            app_rad = np.hypot(a[:, 4] - r0, a[:, 5] - z0)

            # BDD 2a: radial envelope. The app reaches the same outer boundary
            # and never overshoots it; it may extend inward of the reference in
            # the axis volume because it seeds down to s = -0.85 (ref: -0.5).
            if app_rad.min() < -1e-9:
                failures.append(f"BDD2 vol{lvol}: negative radius")
            if app_rad.max() > ref_rad.max() * (1.0 + ENV_OVERSHOOT):
                failures.append(
                    f"BDD2 vol{lvol}: app overshoots outer envelope "
                    f"{app_rad.max():.4f} > {ref_rad.max():.4f}")
            if app_rad.max() < ref_rad.max() * (1.0 - ENV_REACH):
                failures.append(
                    f"BDD2 vol{lvol}: app fails to reach outer envelope "
                    f"{app_rad.max():.4f} < {ref_rad.max():.4f}")

            # BDD 2b: precise geometry on the shared outermost surface.
            outer = int(np.unique(a[:, 0]).max())
            m = a[a[:, 0] == outer]
            nn = median_nearest(np.column_stack([m[:, 4], m[:, 5]]), ref_pts)
            if not nn < NN_MEDIAN_TOL:
                failures.append(
                    f"BDD2 vol{lvol}: outer-surface nn median {nn:.3e} "
                    f">= {NN_MEDIAN_TOL:.1e} m")
            print(f"BDD2 vol{lvol}: outer nn median = {nn:.3e} m, "
                  f"app_rad[{app_rad.min():.4f},{app_rad.max():.4f}] "
                  f"ref_rad[{ref_rad.min():.4f},{ref_rad.max():.4f}]")

            # BDD 3: rotational transform in the 10%-expanded reference band
            # with the same (positive) sign -- catches a reversed orientation.
            iota = load_iota(work, lvol)
            ref_iota = ref[f"iota_vol{i}"]
            w = ref_iota.max() - ref_iota.min()
            lo, hi = ref_iota.min() - IOTA_EXPAND * w, ref_iota.max() + IOTA_EXPAND * w
            if not (iota >= lo).all() or not (iota <= hi).all():
                failures.append(
                    f"BDD3 vol{lvol}: iota {iota} outside [{lo:.5f},{hi:.5f}]")
            if np.sign(iota).min() != np.sign(ref_iota).max():
                failures.append(
                    f"BDD3 vol{lvol}: iota sign {np.sign(iota)} != "
                    f"reference {np.sign(ref_iota)}")
            print(f"BDD3 vol{lvol}: iota [{iota.min():.5f},{iota.max():.5f}] "
                  f"ref [{ref_iota.min():.5f},{ref_iota.max():.5f}]")

        # BDD 4: implicit midpoint is second order, so the final-section shift
        # falls by about four when the step is halved.
        pts = {}
        for nsteps in (256, 512, 1024):
            run_app(binary, work, h5, nturns=3, nsteps=nsteps, nseeds=4,
                    volumes=[2])
            pts[nsteps] = final_section_point(load_volume(work, 2), 2)

        def shift(n1, n2):
            ds = abs(pts[n1][0] - pts[n2][0])
            dth = abs((pts[n1][1] - pts[n2][1] + np.pi) % (2 * np.pi) - np.pi)
            return np.hypot(ds, dth)

        d1, d2 = shift(256, 512), shift(512, 1024)
        if not d1 < 1.0e-3:
            failures.append(f"BDD4: coarse shift {d1:.3e} not small")
        if not d2 < 0.40 * d1:
            failures.append(f"BDD4: shift not converging {d2:.3e} vs {d1:.3e}")
        print(f"BDD4: shift 256->512 = {d1:.3e}, 512->1024 = {d2:.3e}")

        seed_path = os.path.join(work, "seeds.dat")
        explicit = np.array([[1, -0.4, 0.1], [1, 0.2, 0.7]])
        np.savetxt(seed_path, explicit, fmt=["%d", "%.17g", "%.17g"])
        run_app(binary, work, h5, nturns=1, nsteps=64, nseeds=1,
                volumes=[2], seed_file=seed_path)
        initial = load_volume(work, 2)
        initial = initial[initial[:, 1] == 0][:, [2, 3]]
        if not np.allclose(initial, explicit[:, 1:], rtol=0.0, atol=1.0e-14):
            failures.append("BDD5: explicit field-line seeds were replaced")
        print("BDD5: explicit field-line seeds preserved")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  " + f)
        sys.exit(1)
    print("\nAll BDD scenarios passed.")


if __name__ == "__main__":
    main()
