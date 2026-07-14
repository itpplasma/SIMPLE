#!/usr/bin/env python3
"""End-to-end SPECTRE validation suite (#442).

Production-scope validation of SPECTRE guiding-center support on the committed
two-volume tok2vol fixture (integ_coords = 6, integmode = 3, crossing_level = 1).
Markers run at facE_al = 4 (875 keV): the guiding-center-defensible energy on
tok2vol -- at the full 3.5 MeV the interior GC degeneracy breaks the canonical
model for a subset of pitches (see DOC/spectre-interface-crossing.md and
test_spectre_gc_rk45).

Scenarios:
  1. Losses + accounting: every marker is confined or terminated early, no marker
     is left unresolved, and confined_fraction.dat closes the same account that
     times_lost.dat reports (two independently maintained counters must agree).
  2. Level-0 and Level-1 crossing maps: the same fixed-seed ensemble under both
     maps closes particle accounting and reproduces bit-for-bit. Their loss
     fractions are reported separately because the maps encode different sheet
     physics and need not produce the same transport.
  3. Classification: the trapped/passing split matches the analytic mirror
     criterion v_par^2 < 2 mu (Bmax - B). With v = v0 at launch this is
     perp_inv = v_perp^2/B > 1/Bmax, checked independently of the code's trap_par
     sign from perp_inv (times_lost.dat) and the printed surface Bmax.
  4. Guiding-center crossing mu: the Level-1 map holds mu fixed within every
     crossing, so a marker's logged mu drifts only by the symplectic GC flow
     between crossings. That GC-side mu scatter is printed and bounded (it
     quantifies the conserved-mu model's self-consistency). The PHYSICAL mu
     scatter against a full-orbit Boris reference is ROADMAP (see the mu note).
  5. Determinism: a fixed-seed run reproduces times_lost.dat bit-for-bit.

Usage: test_spectre_validation.py <simple.x> <tok2vol.h5>
"""
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

TRACE_TIME = 1.0e-4
NPART = 64
FACE_AL = 4.0
NPOIPER2 = 128

# spectre_crossing_events.dat columns.
C_PART, C_TYPE, C_MU = 0, 3, 10
TYPE_CROSSING, TYPE_REFLECTION = 1, 2

# The Level-1 map conserves mu within each crossing exactly; the residual is the
# symplectic GC flow's between-crossing mu drift (measured ~4e-5 here).
MU_SCATTER_TOL = 5.0e-3
# The mirror criterion is exact, so the analytic split and trap_par must agree.
CLASS_MATCH_MIN = 1.0


def write_input(path, h5, crossing_level=1, ntestpart=NPART,
                trace_time=TRACE_TIME, sbeg=0.5, ntimstep=100,
                npoiper2=NPOIPER2, relerr="1d-11", face_al=FACE_AL,
                ran_seed=12345):
    lines = [
        "&config",
        f"  trace_time = {trace_time}",
        f"  sbeg = {sbeg}",
        f"  ntestpart = {ntestpart}",
        f"  ntimstep = {ntimstep}",
        f"  npoiper2 = {npoiper2}",
        f"  relerr = {relerr}",
        f"  facE_al = {face_al}",
        f"  field_input = '{h5}'",
        "  integ_coords = 6",
        "  integmode = 3",
        "  spectre_ncon_phi = 32",
        f"  crossing_level = {crossing_level}",
        "  deterministic = .True.",
        f"  ran_seed = {ran_seed}",
        "/",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run(binary, workdir, h5, **kwargs):
    write_input(os.path.join(workdir, "simple.in"), h5, **kwargs)
    proc = subprocess.run([binary, "simple.in"], cwd=workdir,
                          capture_output=True, text=True, timeout=280)
    if proc.returncode != 0:
        raise RuntimeError(f"simple.x failed:\n{proc.stdout[-3000:]}\n"
                           f"{proc.stderr[-1500:]}")
    return proc.stdout


def split_confined(times_lost, trace_time):
    confined = times_lost >= trace_time * (1.0 - 1e-9)
    terminated = (times_lost > 0.0) & (times_lost < trace_time * (1.0 - 1e-9))
    unresolved = times_lost <= 0.0
    return confined, terminated, unresolved


def losses_and_accounting(binary, h5, failures):
    """Run the headline L1 ensemble and check the account. Returns the loss
    fraction so the Level-0 comparison reuses this run."""
    with tempfile.TemporaryDirectory() as work:
        run(binary, work, h5, crossing_level=1)
        tl = np.loadtxt(os.path.join(work, "times_lost.dat"))
        cf = np.loadtxt(os.path.join(work, "confined_fraction.dat"))
        ev_path = os.path.join(work, "spectre_crossing_events.dat")
        crossings_written = os.path.exists(ev_path) and os.path.getsize(ev_path) > 0

    if len(tl) != NPART:
        failures.append(f"accounting: times_lost.dat has {len(tl)} rows != {NPART}")
    if not np.all(np.isfinite(tl[:, 1])):
        failures.append("accounting: non-finite loss time (unresolved marker)")

    confined, terminated, unresolved = split_confined(tl[:, 1], TRACE_TIME)
    n_conf, n_term, n_unres = (int(confined.sum()), int(terminated.sum()),
                               int(unresolved.sum()))
    if n_unres != 0:
        failures.append(f"accounting: {n_unres} markers left unresolved "
                        f"(times_lost <= 0)")
    if n_conf + n_term != NPART:
        failures.append(f"accounting: {n_conf} confined + {n_term} terminated "
                        f"!= {NPART} markers")

    # confined_fraction.dat is filled by an independent per-timestep confined
    # counter; its final value must equal the times_lost confined count. Breaking
    # the loss tally (counting a lost marker as confined) makes these disagree.
    frac_final = cf[-1, 1] + cf[-1, 2]
    if abs(frac_final - n_conf / NPART) > 0.5 / NPART + 1e-9:
        failures.append(f"accounting: confined_fraction {frac_final:.4f} != "
                        f"times_lost confined {n_conf / NPART:.4f}")
    if int(cf[-1, 3]) != NPART:
        failures.append(f"accounting: confined_fraction N = {int(cf[-1, 3])}")
    if not crossings_written:
        failures.append("accounting: no spectre_crossing_events.dat produced")

    print(f"accounting: confined={n_conf} terminated={n_term} unresolved="
          f"{n_unres} (sum={n_conf + n_term}={NPART}) "
          f"cf_final={frac_final:.3f}=({n_conf}/{NPART})")
    return n_term / NPART


def crossing_map_accounting(binary, h5, p1, failures):
    with tempfile.TemporaryDirectory() as work:
        run(binary, work, h5, crossing_level=0)
        tl0 = np.loadtxt(os.path.join(work, "times_lost.dat"))
        raw0 = Path(work, "times_lost.dat").read_bytes()
    confined0, terminated0, unresolved0 = split_confined(tl0[:, 1], TRACE_TIME)
    n_conf0, n_term0, n_unres0 = (int(confined0.sum()), int(terminated0.sum()),
                                  int(unresolved0.sum()))
    if n_conf0 + n_term0 != NPART or n_unres0 != 0:
        failures.append(f"crossing maps: Level-0 account is {n_conf0} confined + "
                        f"{n_term0} terminated + {n_unres0} unresolved")

    with tempfile.TemporaryDirectory() as work:
        run(binary, work, h5, crossing_level=0)
        raw0_repeat = Path(work, "times_lost.dat").read_bytes()
    if raw0 != raw0_repeat:
        failures.append("crossing maps: Level-0 times_lost.dat is not "
                        "bit-identical under a fixed seed")

    p0 = n_term0 / NPART
    print(f"crossing maps: loss_fraction L1={p1:.3f} L0={p0:.3f} "
          f"Level-0_account={n_conf0}+{n_term0}+{n_unres0}={NPART} "
          f"Level-0_repeat={'bit-identical' if raw0 == raw0_repeat else 'DIFFERS'}")


def classification(binary, h5, failures):
    with tempfile.TemporaryDirectory() as work:
        out = run(binary, work, h5, ntestpart=128, trace_time=1.0e-6, ntimstep=5,
                  relerr="1d-8")
        tl = np.loadtxt(os.path.join(work, "times_lost.dat"))
        start = np.loadtxt(os.path.join(work, "start.dat"))

    m = re.search(r"bmin = *([0-9.E+-]+) *bmax = *([0-9.E+-]+)", out)
    if not m:
        failures.append("classification: surface bmin/bmax line missing")
        return
    bmin, bmax = float(m.group(1)), float(m.group(2))

    trap_par, perp_inv, lam = tl[:, 2], tl[:, 4], start[:, 4]
    code_trapped = trap_par > 0.0
    mirror_trapped = perp_inv > 1.0 / bmax   # v_par^2 < 2 mu (Bmax - B), v = v0
    match = float(np.mean(code_trapped == mirror_trapped))
    if match < CLASS_MATCH_MIN:
        failures.append(f"classification: mirror criterion matches trap_par on "
                        f"{match:.1%} of markers (< {CLASS_MATCH_MIN:.0%})")
    if not (code_trapped.any() and (~code_trapped).any()):
        failures.append("classification: split degenerate (one class empty)")

    # perp_inv encodes a real on-surface |B| = (1 - lambda^2)/perp_inv in
    # [Bmin, Bmax] (independent physical consistency of the sampled markers).
    b_marker = (1.0 - lam ** 2) / perp_inv
    if not np.all(b_marker <= bmax * (1.0 + 1e-6)):
        failures.append("classification: reconstructed |B| exceeds surface Bmax")
    if not np.all(b_marker >= bmin * (1.0 - 1e-6)):
        failures.append("classification: reconstructed |B| below surface Bmin")

    print(f"classification: trapped={int(code_trapped.sum())}/{len(tl)} "
          f"(frac={code_trapped.mean():.3f}) mirror_match={match:.1%} "
          f"Bmax={bmax:.3e} Bmin={bmin:.3e}")


def gc_mu_conservation(binary, h5, failures):
    # Small dedicated run near the interface: a few markers cross repeatedly,
    # keeping the crossing log small enough to load.
    with tempfile.TemporaryDirectory() as work:
        run(binary, work, h5, ntestpart=8, sbeg=0.9, trace_time=5.0e-6,
            ntimstep=50, npoiper2=256)
        ev = np.loadtxt(os.path.join(work, "spectre_crossing_events.dat"))
    ev = ev.reshape(1, -1) if ev.ndim == 1 else ev
    events = ev[np.isin(ev[:, C_TYPE], (TYPE_CROSSING, TYPE_REFLECTION))]
    if len(events) < 10:
        failures.append(f"mu: only {len(events)} crossing events (< 10)")
        return

    scatter = 0.0
    for p in np.unique(events[:, C_PART]):
        mus = events[events[:, C_PART] == p][:, C_MU]
        ref = float(np.median(mus))
        if ref > 0.0:
            scatter = max(scatter, float(np.max(np.abs(mus - ref)) / ref))
    if not scatter < MU_SCATTER_TOL:
        failures.append(f"mu: GC-map mu scatter {scatter:.3e} >= "
                        f"{MU_SCATTER_TOL:.0e} (conserved-mu model self-broken)")
    print(f"mu: crossing_events={len(events)} GC-side mu_scatter={scatter:.3e} "
          f"(conserved-mu model; the map holds mu fixed per crossing, this is "
          f"the between-crossing symplectic drift)")
    print("mu: ROADMAP -- the PHYSICAL mu scatter needs a full-orbit Boris trace "
          "across a SPECTRE interface; orbit_model=7 evaluates the Boozer/"
          "chartmap flux potential through ref_coords and has no SPECTRE "
          "per-volume Cartesian map, so Boris-on-SPECTRE is not wired.")


def determinism(binary, h5, failures):
    with tempfile.TemporaryDirectory() as w1:
        run(binary, w1, h5, crossing_level=1, trace_time=2.0e-5)
        tl1 = (Path(w1) / "times_lost.dat").read_bytes()
    with tempfile.TemporaryDirectory() as w2:
        run(binary, w2, h5, crossing_level=1, trace_time=2.0e-5)
        tl2 = (Path(w2) / "times_lost.dat").read_bytes()
    identical = tl1 == tl2
    if not identical:
        failures.append("determinism: times_lost.dat not reproducible with a "
                        "fixed seed")
    print(f"determinism: times_lost.dat "
          f"{'bit-identical' if identical else 'DIFFERS'} across fixed-seed runs")


def main():
    binary, tok2vol = sys.argv[1], sys.argv[2]
    failures = []

    p1 = losses_and_accounting(binary, tok2vol, failures)
    crossing_map_accounting(binary, tok2vol, p1, failures)
    classification(binary, tok2vol, failures)
    gc_mu_conservation(binary, tok2vol, failures)
    determinism(binary, tok2vol, failures)

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  " + f)
        sys.exit(1)
    print("\nAll SPECTRE validation scenarios passed.")


if __name__ == "__main__":
    main()
