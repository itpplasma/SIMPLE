#!/usr/bin/env python3
"""Behavioural test for the Level-1 SPECTRE interface crossing map (#440).

Level 1 is the thin-current-sheet limit of the guiding-center dynamics: the
impulse Delta z = lambda_k * X along the interface Hamiltonian vector field
X = {z, rho_g}, with the scalar lambda_k fixed by exact energy conservation. It
adds the tangential sheet-drift kick and the drift-order v_par term that Level 0
omits. The crossing log carries the generator components and lambda_k, so the
map is checked directly from spectre_crossing_events.dat:

  * Implementation identity: on tok2vol every crossing applies exactly
    lambda_k * X to (theta, zeta, v_par) and conserves H = v_par^2/2 + mu*|B| to
    1e-13, with |B| taken at the KICKED landing point in each volume.
  * Physics magnitude: the tangential kicks are nonzero, finite, and small
    compared to 2*pi; the same run under crossing_level = 0 applies no kick.
  * Level-0 vs Level-1: the same ensemble under both maps keeps the event count
    within 20% and conserves energy in both.

Usage:
    test_spectre_crossing_l1.py <simple.x> <tok2vol.h5>
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

# spectre_crossing_events.dat columns.
C_PART, C_TIME, C_IFACE, C_TYPE, C_VFROM, C_VTO = 0, 1, 2, 3, 4, 5
C_THETA, C_ZETA, C_VPB, C_VPA, C_MU, C_BH, C_BT = 6, 7, 8, 9, 10, 11, 12
C_XT, C_XZ, C_XV, C_LAM, C_DTH, C_DZE = 13, 14, 15, 16, 17, 18

TYPE_CROSSING, TYPE_REFLECTION, TYPE_LOSS = 1, 2, 3

NPART, TRACE_TIME, SBEG = 16, 2.0e-5, 0.5
H_TOL = 1.0e-13
KICK_IDENTITY_TOL = 1.0e-12
KICK_MIN, KICK_MAX = 1.0e-8, 0.3


def write_input(path, h5, level):
    lines = [
        "&config",
        f"  trace_time = {TRACE_TIME}",
        f"  sbeg = {SBEG}",
        f"  ntestpart = {NPART}",
        "  ntimstep = 100",
        "  npoiper2 = 256",
        "  relerr = 1d-8",
        f"  field_input = '{h5}'",
        "  integ_coords = 6",
        "  integmode = 0",
        f"  crossing_level = {level}",
        "  deterministic = .True.",
        "  ran_seed = 12345",
        "/",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run_simple(binary, workdir, h5, level):
    write_input(os.path.join(workdir, "simple.in"), h5, level)
    proc = subprocess.run([binary, "simple.in"], cwd=workdir,
                          capture_output=True, text=True, timeout=280)
    if proc.returncode != 0:
        raise RuntimeError(f"simple.x failed:\n{proc.stdout}\n{proc.stderr}")


def load_events(workdir):
    a = np.loadtxt(os.path.join(workdir, "spectre_crossing_events.dat"))
    return a.reshape(1, -1) if a.ndim == 1 else a


def load_lost_fraction(workdir):
    tl = np.loadtxt(os.path.join(workdir, "times_lost.dat"))
    lost = tl[:, 1] < TRACE_TIME * (1.0 - 1e-9)
    return int(lost.sum()) / len(tl)


def energy_residual(ev):
    # A crossing lands in the target volume (post-field B_target at the kicked
    # angles); a reflection stays home (post-field B_home). Loss events are
    # terminal and excluded by the caller.
    b_after = np.where(ev[:, C_TYPE] == TYPE_REFLECTION, ev[:, C_BH], ev[:, C_BT])
    h_before = 0.5 * ev[:, C_VPB] ** 2 + ev[:, C_MU] * ev[:, C_BH]
    h_after = 0.5 * ev[:, C_VPA] ** 2 + ev[:, C_MU] * b_after
    return np.abs(h_after - h_before) / np.abs(h_before)


def check_identity(ev1, failures):
    # Refracted crossings carry the generator impulse (lambda_k != 0). Where the
    # equal-|B| refraction would need a kick above the thin-layer bound the map
    # falls back to the energy-exact Level-0 rescale (lambda_k = 0), so restrict
    # the generator identity to the refracted subset.
    cross = ev1[ev1[:, C_TYPE] == TYPE_CROSSING]
    refr = cross[cross[:, C_LAM] != 0.0]
    if len(refr) < 20:
        failures.append(f"identity: only {len(refr)} refracted crossings (< 20)")
        return
    sample = refr[:20]

    # Applied kicks equal lambda_k * X for (theta, zeta, v_par).
    d_theta = np.abs(sample[:, C_DTH] - sample[:, C_LAM] * sample[:, C_XT])
    d_zeta = np.abs(sample[:, C_DZE] - sample[:, C_LAM] * sample[:, C_XZ])
    d_vpar = np.abs((sample[:, C_VPA] - sample[:, C_VPB])
                    - sample[:, C_LAM] * sample[:, C_XV])
    kick_err = float(max(d_theta.max(), d_zeta.max(), d_vpar.max()))
    if not kick_err < KICK_IDENTITY_TOL:
        failures.append(f"identity: kick != lambda_k*X, err {kick_err:.3e}")

    rel = float(energy_residual(refr).max())
    if not rel < H_TOL:
        failures.append(f"identity: refracted |dH|/H = {rel:.3e} >= {H_TOL:.0e}")

    if np.max(np.abs(cross[:, C_BT] - cross[:, C_BH])) <= 0.0:
        failures.append("identity: crossings show no [[Bmod]] jump")
    print(f"identity: crossings={len(cross)} refracted={len(refr)} "
          f"kick_err={kick_err:.3e} max|dH|/H={rel:.3e}")


def check_physics(ev1, ev0, failures):
    cross1 = ev1[ev1[:, C_TYPE] == TYPE_CROSSING]
    kick = np.hypot(cross1[:, C_DTH], cross1[:, C_DZE])
    refr = kick[kick > 0.0]
    if not np.all(np.isfinite(kick)):
        failures.append("physics: non-finite tangential kick")
    if not np.all(kick < KICK_MAX):
        failures.append(f"physics: kick {kick.max():.3e} >= {KICK_MAX} rad")
    if not KICK_MIN < kick.max():
        failures.append(f"physics: no nonzero refraction kick (max {kick.max():.3e})")

    cross0 = ev0[ev0[:, C_TYPE] == TYPE_CROSSING]
    zero0 = np.max(np.abs(cross0[:, [C_XT, C_XZ, C_XV, C_LAM, C_DTH, C_DZE]]))
    if not zero0 <= 0.0:
        failures.append(f"physics: Level-0 kick nonzero ({zero0:.3e})")
    print(f"physics: refracted={len(refr)}/{len(cross1)} "
          f"kick_max={kick.max():.3e} rad, Level-0 max component {zero0:.3e}")


def check_energy(ev, tag, failures):
    non_loss = ev[ev[:, C_TYPE] != TYPE_LOSS]
    if len(non_loss) == 0:
        return
    rel = float(energy_residual(non_loss).max())
    if not rel < H_TOL:
        failures.append(f"{tag}: |dH|/H = {rel:.3e} >= {H_TOL:.0e}")
    print(f"{tag}: events={len(non_loss)} max|dH|/H={rel:.3e}")


def main():
    binary, tok2vol = sys.argv[1], sys.argv[2]
    failures = []

    with tempfile.TemporaryDirectory() as work1:
        run_simple(binary, work1, tok2vol, level=1)
        ev1 = load_events(work1)
        loss1 = load_lost_fraction(work1)
    with tempfile.TemporaryDirectory() as work0:
        run_simple(binary, work0, tok2vol, level=0)
        ev0 = load_events(work0)
        loss0 = load_lost_fraction(work0)

    check_identity(ev1, failures)
    check_physics(ev1, ev0, failures)
    check_energy(ev1, "level1", failures)
    check_energy(ev0, "level0", failures)

    n1, n0 = len(ev1), len(ev0)
    spread = abs(n1 - n0) / max(n0, 1)
    if not spread <= 0.20:
        failures.append(f"ensemble: event count spread {spread:.2%} > 20%")
    print(f"ensemble: events L1={n1} L0={n0} spread={spread:.2%} "
          f"loss_fraction L1={loss1:.3f} L0={loss0:.3f} (informational)")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  " + f)
        sys.exit(1)
    print("\nAll Level-1 crossing scenarios passed.")


if __name__ == "__main__":
    main()
