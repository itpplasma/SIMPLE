#!/usr/bin/env python3
"""Behavioural test for symplectic multi-volume interface crossing (#441).

Runs simple.x on the committed tok2vol fixture (integ_coords=6, integmode=3)
and checks the exact-landing crossing pipeline end to end:

  * Landing accuracy: hundreds of symplectic crossings land on the interface
    with |rho_g - k| < 1e-8 (read from the sympl_landing_stats stdout line);
    zero CROSS_STOP events in every run.
  * Reflection: deeply-trapped markers mirror at the interface (events logged),
    stay in their volume, conserve H = v_par^2/2 + mu*B exactly across the
    event, and the orbit continues afterwards.
  * Cross-path consistency: the same 32-marker ensemble under integmode 3 and
    integmode 0 (RK45 + crossing) loses a fraction of markers that agrees
    within Monte-Carlo error.
  * Step-halving: a crossing-rich single orbit re-run with doubled npoiper2
    converges at the scheme's order (delta ratio reported and bounded).

Usage:
    test_spectre_sympl_crossing.py <simple.x> <tok2vol.h5>
"""
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
from netCDF4 import Dataset

# spectre_crossing_events.dat columns.
C_PART, C_TIME, C_IFACE, C_TYPE, C_VFROM, C_VTO = 0, 1, 2, 3, 4, 5
C_THETA, C_ZETA, C_VPB, C_VPA, C_MU, C_BH, C_BT = 6, 7, 8, 9, 10, 11, 12

TYPE_CROSSING, TYPE_REFLECTION, TYPE_LOSS, TYPE_STOP = 1, 2, 3, 4

LANDING_TOL = 1.0e-8
H_REFL_TOL = 1.0e-12
MIN_LANDINGS = 200
HALVING_RATIO_LO, HALVING_RATIO_HI = 3.0, 6.0


def write_input(path, h5, integmode, npart, trace_time, sbeg, npoiper2,
                relerr, face_al=1.0, ntimstep=100):
    lines = [
        "&config",
        f"  trace_time = {trace_time}",
        f"  sbeg = {sbeg}",
        f"  ntestpart = {npart}",
        f"  ntimstep = {ntimstep}",
        f"  npoiper2 = {npoiper2}",
        f"  relerr = {relerr}",
        f"  facE_al = {face_al}",
        f"  field_input = '{h5}'",
        "  integ_coords = 6",
        f"  integmode = {integmode}",
        "  output_orbits_macrostep = .True.",
        "  deterministic = .True.",
        "  ran_seed = 12345",
        "/",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run_simple(binary, workdir, **kwargs):
    write_input(os.path.join(workdir, "simple.in"), **kwargs)
    proc = subprocess.run([binary, "simple.in"], cwd=workdir,
                          capture_output=True, text=True, timeout=280)
    if proc.returncode != 0:
        raise RuntimeError(f"simple.x failed:\n{proc.stdout[-3000:]}\n"
                           f"{proc.stderr[-2000:]}")
    return proc.stdout


def load_events(workdir):
    a = np.loadtxt(os.path.join(workdir, "spectre_crossing_events.dat"))
    if a.size == 0:
        return np.empty((0, 19))
    return a.reshape(1, -1) if a.ndim == 1 else a


def load_times_lost(workdir):
    return np.loadtxt(os.path.join(workdir, "times_lost.dat"))


def parse_landing_stats(stdout):
    m = re.search(r"sympl_landing_stats: count=\s*(\d+)\s+max_resid=\s*"
                  r"([0-9.E+-]+)\s+cross_stop=\s*(\d+)", stdout)
    if not m:
        raise RuntimeError("sympl_landing_stats line missing from stdout")
    return int(m.group(1)), float(m.group(2)), int(m.group(3))


def check_stops(ev, stops, tag, failures):
    n_stop_log = int((ev[:, C_TYPE] == TYPE_STOP).sum()) if len(ev) else 0
    if stops != 0 or n_stop_log != 0:
        failures.append(f"{tag}: CROSS_STOP events ({stops}, log {n_stop_log})")


def check_landing_and_reflection(binary, h5, failures):
    # 875 keV alphas: guiding-center-valid on tok2vol. At the full 3.5 MeV the
    # interior GC degeneracy (Bstar_par = 0 mid-volume) breaks the canonical
    # model for a subset of pitches and those orbits stop, correctly.
    trace_time, npart = 2.0e-5, 16
    with tempfile.TemporaryDirectory() as work:
        out = run_simple(binary, work, h5=h5, integmode=3, npart=npart,
                         trace_time=trace_time, sbeg=0.5, npoiper2=256,
                         relerr="1d-12", face_al=4.0)
        ev = load_events(work)
        tl = load_times_lost(work)
    landings, max_resid, stops = parse_landing_stats(out)
    check_stops(ev, stops, "landing", failures)

    if landings < MIN_LANDINGS:
        failures.append(f"landing: only {landings} landings (< {MIN_LANDINGS})")
    if not max_resid < LANDING_TOL:
        failures.append(f"landing: max |rho_g - k| = {max_resid:.3e} "
                        f">= {LANDING_TOL:.0e}")
    print(f"landing: landings={landings} max|rho_g - k|={max_resid:.3e}")

    refl = ev[ev[:, C_TYPE] == TYPE_REFLECTION]
    if len(refl) < 1:
        failures.append("reflection: no reflection events under integmode 3")
        return
    if not np.all(refl[:, C_VFROM] == refl[:, C_VTO]):
        failures.append("reflection: marker switched volume")
    h_before = 0.5*refl[:, C_VPB]**2 + refl[:, C_MU]*refl[:, C_BH]
    h_after = 0.5*refl[:, C_VPA]**2 + refl[:, C_MU]*refl[:, C_BH]
    rel = float(np.max(np.abs(h_after - h_before)/np.abs(h_before)))
    if not rel < H_REFL_TOL:
        failures.append(f"reflection: |dH|/H = {rel:.3e} >= {H_REFL_TOL:.0e}")

    # The orbit continues after reflecting: some reflecting marker survives
    # beyond its last reflection (confined for the full trace or lost later).
    continues = False
    for p in np.unique(refl[:, C_PART]).astype(int):
        t_last = refl[refl[:, C_PART] == p][:, C_TIME].max()
        t_end = float(tl[tl[:, 0] == p][0, 1])
        if t_end >= trace_time*(1.0 - 1e-9) or t_end > 2.0*t_last:
            continues = True
    if not continues:
        failures.append("reflection: no reflecting marker continued its orbit")
    print(f"reflection: events={len(refl)} max|dH|/H={rel:.3e} "
          f"continues={continues}")


def loss_fraction(workdir, trace_time):
    tl = load_times_lost(workdir)
    lost = tl[:, 1] < trace_time*(1.0 - 1e-9)
    return float(lost.sum())/len(tl), len(tl)


def check_cross_path(binary, h5, failures):
    # 875 keV: the most guiding-center-defensible energy on tok2vol at which
    # the RK45 reference path still integrates (its odeint aborts at lower
    # energies); at full 3.5 MeV the two GC formulations differ at O(rho*) ~ 1
    # and their loss statistics are not comparable.
    trace_time, npart = 2.0e-5, 32
    with tempfile.TemporaryDirectory() as work:
        out = run_simple(binary, work, h5=h5, integmode=3, npart=npart,
                         trace_time=trace_time, sbeg=0.5, npoiper2=256,
                         relerr="1d-12", face_al=4.0)
        ev3 = load_events(work)
        _, _, stops = parse_landing_stats(out)
        check_stops(ev3, stops, "cross-path sympl", failures)
        p3, n = loss_fraction(work, trace_time)
    with tempfile.TemporaryDirectory() as work:
        run_simple(binary, work, h5=h5, integmode=0, npart=npart,
                   trace_time=trace_time, sbeg=0.5, npoiper2=256,
                   relerr="1d-8", face_al=4.0)
        p0, _ = loss_fraction(work, trace_time)

    p_mean = 0.5*(p3 + p0)
    sigma = np.sqrt(max(p_mean*(1.0 - p_mean), 1.0/n)/n)
    if not abs(p3 - p0) < 3.0*sigma:
        failures.append(f"cross-path: loss fractions differ beyond MC error: "
                        f"sympl {p3:.3f} vs RK45 {p0:.3f} (3 sigma = "
                        f"{3.0*sigma:.3f})")
    print(f"cross-path: loss fraction sympl={p3:.3f} RK45={p0:.3f} "
          f"3sigma={3.0*sigma:.3f} (N={n})")


def orbit_states(workdir):
    ds = Dataset(os.path.join(workdir, "orbits.nc"))
    z = np.stack([np.array(ds[v][:]) for v in
                  ("s", "theta", "phi", "p_abs", "v_par")])
    ds.close()
    return z  # (component, timestep, particle)


def check_step_halving(binary, h5, failures):
    trace_time, npart = 5.0e-5, 8
    states = []
    crossers = []
    confined = []
    for npoiper2 in (1024, 2048, 4096):
        with tempfile.TemporaryDirectory() as work:
            out = run_simple(binary, work, h5=h5, integmode=3, npart=npart,
                             trace_time=trace_time, sbeg=0.97,
                             npoiper2=npoiper2, relerr="1d-12", face_al=50.0)
            ev = load_events(work)
            _, _, stops = parse_landing_stats(out)
            check_stops(ev, stops, f"halving np={npoiper2}", failures)
            states.append(orbit_states(work))
            cross = ev[ev[:, C_TYPE] == TYPE_CROSSING] if len(ev) else ev
            crossers.append(set(cross[:, C_PART].astype(int)) if len(cross)
                            else set())
            tl = load_times_lost(work)
            confined.append(set(tl[:, 0][tl[:, 1] >= trace_time*(1 - 1e-9)]
                                .astype(int)))

    usable = (crossers[0] & crossers[1] & crossers[2]
              & confined[0] & confined[1] & confined[2])
    if not usable:
        failures.append(f"halving: no marker crossing and confined in all "
                        f"runs (crossers {sorted(crossers[0])}, confined "
                        f"{sorted(confined[0])})")
        return
    ratios = []
    for p in sorted(usable):
        z0, z1, z2 = (s[:, -1, p - 1] for s in states)
        d1 = np.linalg.norm(z0 - z1)
        d2 = np.linalg.norm(z1 - z2)
        ratios.append((p, d1, d2, d1/d2 if d2 > 0 else np.inf))
    # The convergence order is read from the marker with the cleanest signal
    # (largest d2 above round-off).
    p, d1, d2, ratio = max(ratios, key=lambda r: r[2])
    if not HALVING_RATIO_LO <= ratio <= HALVING_RATIO_HI:
        failures.append(f"halving: convergence ratio {ratio:.2f} outside "
                        f"[{HALVING_RATIO_LO}, {HALVING_RATIO_HI}] "
                        f"(marker {p}, d1={d1:.3e}, d2={d2:.3e})")
    print(f"halving: markers={sorted(usable)} marker {p}: "
          f"d(h,h/2)={d1:.3e} d(h/2,h/4)={d2:.3e} ratio={ratio:.2f}")


def main():
    binary, tok2vol = sys.argv[1], sys.argv[2]
    failures = []

    check_landing_and_reflection(binary, tok2vol, failures)
    check_cross_path(binary, tok2vol, failures)
    check_step_halving(binary, tok2vol, failures)

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  " + f)
        sys.exit(1)
    print("\nAll symplectic crossing scenarios passed.")


if __name__ == "__main__":
    main()
