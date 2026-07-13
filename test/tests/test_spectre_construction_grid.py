#!/usr/bin/env python3
"""Behavioural test for the SPECTRE construction-grid memory knob (#442).

The per-volume Meiss construction (integmode > 0) allocates rank-3 arrays plus
quintic batch splines of size Mvol*n_r*n_th*n_phi -- ~1.8 GB at the historical
48x48x32 on the two-volume tok2vol fixture. This test covers the config knob
(spectre_ncon_r/th/phi) and the axisymmetric phi auto-clamp:

  * Clamp triggers: tok2vol has no toroidal harmonics, so the default
    (spectre_ncon_phi = -1) auto-clamps the phi grid, printed as
    'axisym=T n_phi=<small> phi_auto=T'.
  * Memory red proof: peak RSS with the clamp is far below the unclamped
    phi=32 build. Forcing spectre_ncon_phi = 32 bypasses the clamp, and RSS
    jumps back by the phi ratio; a coarse r,th grid drops it further still.
  * Construction identity: a coarse-grid run still reproduces the per-volume
    canonical chart -- confined low-energy markers conserve the guiding-center
    energy p^2 to a looser bound (a broken construction drifts O(0.1)).

Usage: test_spectre_construction_grid.py <simple.x> <tok2vol.h5>
"""
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

# Construction identity: confined-marker p^2 drift. The coarse-grid Meiss chart
# stays at ~2.5e-10 (measured); a broken construction drifts to O(0.1).
IDENTITY_TOL = 1.0e-6
# Clamp must cut phi well below the historical 32 to save the phi memory factor.
CLAMP_PHI_MAX = 16
# The clamp must cut peak RSS to a clear fraction of the unclamped phi=32 build.
CLAMP_RSS_FRACTION = 0.6


def write_input(path, h5, integmode, ntestpart, trace_time, ntimstep, npoiper2,
                relerr, face_al, extra=()):
    lines = [
        "&config",
        f"  trace_time = {trace_time}",
        "  sbeg = 0.5",
        f"  ntestpart = {ntestpart}",
        f"  ntimstep = {ntimstep}",
        f"  npoiper2 = {npoiper2}",
        f"  relerr = {relerr}",
        f"  facE_al = {face_al}",
        f"  field_input = '{h5}'",
        "  integ_coords = 6",
        f"  integmode = {integmode}",
        "  deterministic = .True.",
        "  ran_seed = 12345",
        *[f"  {line}" for line in extra],
        "/",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run_peak_rss(binary, workdir):
    """Run simple.x and return (peak_rss_MB, stdout). VmHWM is the kernel's
    per-process resident high-water mark, so polling it captures the transient
    construction peak regardless of read timing."""
    proc = subprocess.Popen([binary, "simple.in"], cwd=workdir,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True)
    peak_kb = 0
    status = f"/proc/{proc.pid}/status"
    while proc.poll() is None:
        try:
            with open(status) as f:
                for line in f:
                    if line.startswith("VmHWM:"):
                        peak_kb = max(peak_kb, int(line.split()[1]))
                        break
        except OSError:
            break
        time.sleep(0.02)
    out = proc.stdout.read()
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"simple.x failed:\n{out[-3000:]}")
    return peak_kb / 1024.0, out


def run_config(binary, h5, integmode=3, ntestpart=2, trace_time=1.0e-6,
               ntimstep=5, npoiper2=256, relerr="1d-8", face_al=1.0, extra=()):
    with tempfile.TemporaryDirectory() as work:
        write_input(os.path.join(work, "simple.in"), h5, integmode, ntestpart,
                    trace_time, ntimstep, npoiper2, relerr, face_al, extra)
        rss, out = run_peak_rss(binary, work)
        tl = os.path.join(work, "times_lost.dat")
        data = np.loadtxt(tl) if os.path.exists(tl) else None
    return rss, out, data


def parse_grid(out):
    m = re.search(r"spectre_construction_grid: n_r=(\d+) n_th=(\d+) "
                  r"axisym=([TF]) n_phi=(\d+) phi_auto=([TF])", out)
    if not m:
        raise RuntimeError("spectre_construction_grid line missing from stdout")
    return dict(n_r=int(m.group(1)), n_th=int(m.group(2)),
                axisym=m.group(3) == "T", n_phi=int(m.group(4)),
                phi_auto=m.group(5) == "T")


def main():
    binary, tok2vol = sys.argv[1], sys.argv[2]
    failures = []

    # Default: axisymmetric tok2vol auto-clamps the phi grid.
    rss_clamped, out_default, _ = run_config(binary, tok2vol)
    grid = parse_grid(out_default)
    if not grid["axisym"]:
        failures.append("clamp: tok2vol not detected as axisymmetric")
    if not grid["phi_auto"]:
        failures.append("clamp: phi_auto false on default config")
    if not grid["n_phi"] <= CLAMP_PHI_MAX:
        failures.append(f"clamp: n_phi {grid['n_phi']} not clamped "
                        f"(> {CLAMP_PHI_MAX})")
    print(f"clamp: n_r={grid['n_r']} n_th={grid['n_th']} axisym={grid['axisym']} "
          f"n_phi={grid['n_phi']} phi_auto={grid['phi_auto']}")

    # Red proof: forcing spectre_ncon_phi = 32 bypasses the clamp and restores
    # the historical full-resolution build; peak RSS jumps by ~the phi ratio.
    rss_full, out_full, _ = run_config(binary, tok2vol,
                                       extra=("spectre_ncon_phi = 32",))
    grid_full = parse_grid(out_full)
    if grid_full["phi_auto"] or grid_full["n_phi"] != 32:
        failures.append(f"red proof: explicit phi=32 not honored "
                        f"(n_phi={grid_full['n_phi']} auto={grid_full['phi_auto']})")
    if not rss_clamped < CLAMP_RSS_FRACTION * rss_full:
        failures.append(f"red proof: clamped RSS {rss_clamped:.0f} MB not below "
                        f"{CLAMP_RSS_FRACTION}x unclamped {rss_full:.0f} MB")
    print(f"memory: clamped(phi={grid['n_phi']})={rss_clamped:.0f} MB  "
          f"unclamped(phi=32)={rss_full:.0f} MB  "
          f"ratio={rss_clamped / rss_full:.2f}")

    # Coarse r,th knob drops peak RSS further below the default.
    rss_coarse, out_coarse, _ = run_config(
        binary, tok2vol,
        extra=("spectre_ncon_r = 16", "spectre_ncon_th = 16"))
    if not rss_coarse < rss_clamped:
        failures.append(f"knob: coarse RSS {rss_coarse:.0f} MB not below "
                        f"default {rss_clamped:.0f} MB")
    print(f"memory: coarse(16x16, phi={parse_grid(out_coarse)['n_phi']})="
          f"{rss_coarse:.0f} MB (default {rss_clamped:.0f} MB)")

    # Construction identity: a coarse-grid chart still conserves the
    # guiding-center energy p^2 for confined low-energy markers.
    _, out_id, data = run_config(
        binary, tok2vol, ntestpart=16, trace_time=1.0e-5, ntimstep=50,
        npoiper2=1024, relerr="1d-12", face_al=50.0,
        extra=("spectre_ncon_r = 16", "spectre_ncon_th = 16"))
    confined = data[:, 1] >= 1.0e-5 * (1.0 - 1e-9)
    if confined.sum() < 8:
        failures.append(f"identity: only {int(confined.sum())} confined markers")
    else:
        drift = float(np.max(np.abs(data[confined, 8] ** 2 - 1.0)))
        if not drift < IDENTITY_TOL:
            failures.append(f"identity: coarse-grid p^2 drift {drift:.3e} "
                            f">= {IDENTITY_TOL:.0e}")
        print(f"identity: coarse-grid confined={int(confined.sum())} "
              f"max|p^2-1|={drift:.3e}")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  " + f)
        sys.exit(1)
    print("\nAll construction-grid scenarios passed.")


if __name__ == "__main__":
    main()
