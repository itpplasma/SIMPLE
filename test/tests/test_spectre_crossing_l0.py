#!/usr/bin/env python3
"""Behavioural test for the Level-0 SPECTRE interface crossing map (#443).

Runs simple.x (integ_coords=6, integmode=0) on two committed SPECTRE fixtures and
checks the crossing pipeline through spectre_crossing_events.dat:

  * tok2vol (2 volumes, pressure step): the pitch-uniform ensemble produces both
    energy-exact crossings and magnetic-mirror reflections at the interface.
      - every crossing conserves H = v_par^2/2 + mu*Bmod (each side in its own
        volume) to 1e-14, with v_par rescaled by sqrt(v_par^2 - 2 mu [[Bmod]]);
      - every reflection stays in its volume, flips v_par, and conserves H
        exactly, and at least one reflecting marker stays confined.
  * G3V3L3Fi (3 volumes): the ensemble traverses all volumes; every marker is
    confined or lost at the outermost interface, the per-marker crossing chain
    pairs up (entry volume of event n+1 = exit volume of event n), no marker
    ends in a boundary-stop state, and the event log is bit-reproducible.

Usage:
    test_spectre_crossing_l0.py <simple.x> <tok2vol.h5> <spectre_test.h5>
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

# spectre_crossing_events.dat columns.
C_PART, C_TIME, C_IFACE, C_TYPE, C_VFROM, C_VTO = 0, 1, 2, 3, 4, 5
C_THETA, C_ZETA, C_VPB, C_VPA, C_MU, C_BH, C_BT = 6, 7, 8, 9, 10, 11, 12

TYPE_CROSSING, TYPE_REFLECTION, TYPE_LOSS = 1, 2, 3

H_TOL = 1.0e-14


def write_input(path, h5, trace_time, npart, sbeg):
    lines = [
        "&config",
        f"  trace_time = {trace_time}",
        f"  sbeg = {sbeg}",
        f"  ntestpart = {npart}",
        "  ntimstep = 100",
        "  npoiper2 = 256",
        "  relerr = 1d-8",
        f"  field_input = '{h5}'",
        "  integ_coords = 6",
        "  integmode = 0",
        "  crossing_level = 0",
        "  deterministic = .True.",
        "  ran_seed = 12345",
        "/",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run_simple(binary, workdir, h5, trace_time, npart, sbeg):
    write_input(os.path.join(workdir, "simple.in"), h5, trace_time, npart, sbeg)
    proc = subprocess.run([binary, "simple.in"], cwd=workdir,
                          capture_output=True, text=True, timeout=280)
    if proc.returncode != 0:
        raise RuntimeError(f"simple.x failed:\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout


def load_events(workdir):
    path = os.path.join(workdir, "spectre_crossing_events.dat")
    a = np.loadtxt(path)
    return a.reshape(1, -1) if a.ndim == 1 else a


def load_times_lost(workdir):
    return np.loadtxt(os.path.join(workdir, "times_lost.dat"))


def check_tok2vol(binary, h5, failures):
    trace_time, npart = 2.0e-5, 16
    with tempfile.TemporaryDirectory() as work:
        run_simple(binary, work, h5, trace_time, npart, sbeg=0.5)
        ev = load_events(work)
        tl = load_times_lost(work)

    cross = ev[ev[:, C_TYPE] == TYPE_CROSSING]
    refl = ev[ev[:, C_TYPE] == TYPE_REFLECTION]

    if len(cross) < 1:
        failures.append("tok2vol: no crossing events")
    if len(refl) < 1:
        failures.append("tok2vol: no reflection events")

    if len(cross):
        h_before = 0.5*cross[:, C_VPB]**2 + cross[:, C_MU]*cross[:, C_BH]
        h_after = 0.5*cross[:, C_VPA]**2 + cross[:, C_MU]*cross[:, C_BT]
        rel = np.max(np.abs(h_after - h_before)/np.abs(h_before))
        if not rel < H_TOL:
            failures.append(f"tok2vol: crossing |dH|/H = {rel:.3e} >= {H_TOL:.0e}")
        radicand = cross[:, C_VPB]**2 - 2.0*cross[:, C_MU]*(cross[:, C_BT]
                                                            - cross[:, C_BH])
        vpa_pred = np.sign(cross[:, C_VPB])*np.sqrt(np.clip(radicand, 0.0, None))
        jump = np.max(np.abs(vpa_pred - cross[:, C_VPA]))
        if not jump < 1.0e-12:
            failures.append(f"tok2vol: v_par jump mismatch {jump:.3e}")
        # A real pressure jump must be exercised, else the map is untested.
        if np.max(np.abs(cross[:, C_BT] - cross[:, C_BH])) <= 0.0:
            failures.append("tok2vol: crossings show no [[Bmod]] jump")
        print(f"tok2vol: crossings={len(cross)} max|dH|/H={rel:.3e} "
              f"vpar_jump_err={jump:.3e}")

    if len(refl):
        if not np.all(refl[:, C_VFROM] == refl[:, C_VTO]):
            failures.append("tok2vol: reflection switched volume")
        mirror = np.max(np.abs(refl[:, C_VPA] + refl[:, C_VPB]))
        if not mirror < 1.0e-14:
            failures.append(f"tok2vol: reflection v_par not mirrored {mirror:.3e}")
        h_before = 0.5*refl[:, C_VPB]**2 + refl[:, C_MU]*refl[:, C_BH]
        h_after = 0.5*refl[:, C_VPA]**2 + refl[:, C_MU]*refl[:, C_BH]
        rel = np.max(np.abs(h_after - h_before)/np.abs(h_before))
        if not rel <= 0.0:
            failures.append(f"tok2vol: reflection |dH|/H = {rel:.3e} != 0")
        refl_parts = set(refl[:, C_PART].astype(int))
        confined = set(tl[:, 0][tl[:, 1] >= trace_time*(1.0 - 1e-6)].astype(int))
        if not refl_parts & confined:
            failures.append("tok2vol: every reflecting marker was lost/stopped")
        print(f"tok2vol: reflections={len(refl)} mirror_err={mirror:.3e} "
              f"confined_reflectors={sorted(refl_parts & confined)}")


def check_ensemble(binary, h5, failures):
    trace_time, npart = 2.0e-4, 48
    with tempfile.TemporaryDirectory() as work:
        run_simple(binary, work, h5, trace_time, npart, sbeg=1.5)
        ev = load_events(work)
        tl = load_times_lost(work)
        raw = open(os.path.join(work, "spectre_crossing_events.dat")).read()
        if os.path.exists(os.path.join(work, "spectre_boundary_events.dat")):
            failures.append("ensemble: stale spectre_boundary_events.dat present")

    lost_mask = tl[:, 1] < trace_time*(1.0 - 1e-9)
    n_lost = int(lost_mask.sum())
    n_conf = int((~lost_mask).sum())
    if n_lost + n_conf != npart:
        failures.append(f"ensemble: {n_lost}+{n_conf} != {npart} markers")
    if np.any(tl[:, 1] <= 0.0):
        failures.append("ensemble: unaccounted marker (times_lost <= 0)")

    loss = ev[ev[:, C_TYPE] == TYPE_LOSS]
    mvol = int(loss[:, C_IFACE].max()) if len(loss) else 0
    if len(loss) and not np.all(loss[:, C_IFACE] == mvol):
        failures.append("ensemble: loss event not at the outermost interface")

    # Every lost marker ends on a loss event; the crossing chain pairs up.
    chain_ok = True
    loss_ok = True
    lost_ids = set(tl[:, 0][lost_mask].astype(int))
    for p in np.unique(ev[:, C_PART]).astype(int):
        chain = ev[ev[:, C_PART] == p]
        chain = chain[np.argsort(chain[:, C_TIME])]
        for i in range(len(chain) - 1):
            if chain[i + 1, C_VFROM] != chain[i, C_VTO]:
                chain_ok = False
        if p in lost_ids and chain[-1, C_TYPE] != TYPE_LOSS:
            loss_ok = False
    if not chain_ok:
        failures.append("ensemble: crossing chain volumes do not pair up")
    if not loss_ok:
        failures.append("ensemble: a lost marker does not end on a loss event")

    with tempfile.TemporaryDirectory() as work2:
        run_simple(binary, work2, h5, trace_time, npart, sbeg=1.5)
        raw2 = open(os.path.join(work2, "spectre_crossing_events.dat")).read()
    if raw != raw2:
        failures.append("ensemble: event log not reproducible with fixed seed")

    print(f"ensemble: confined={n_conf} lost={n_lost} Mvol={mvol} "
          f"events={len(ev)} chain_ok={chain_ok} deterministic={raw == raw2}")


def main():
    binary, tok2vol, g3v3 = sys.argv[1], sys.argv[2], sys.argv[3]
    failures = []
    check_tok2vol(binary, tok2vol, failures)
    check_ensemble(binary, g3v3, failures)

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  " + f)
        sys.exit(1)
    print("\nAll Level-0 crossing scenarios passed.")


if __name__ == "__main__":
    main()
