#!/usr/bin/env python3
"""Compare SPECTRE GC/hybrid/FO cost with Boozer GC on the source VMEC field.

`run NAME INTEGMODE NPOIPER2 MARKERS` traces one Boozer GC case through the
build's simple.x. `report SPECTRE_DIR BOOZER_DIR OUTPUT` normalizes all cases
to core-seconds per simulated particle-second and writes summary.json.
SPECTRE cases come from run.py with matched starts.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

WOUT_URL = ("https://github.com/hiddenSymmetries/simsopt/raw/master/tests/"
            "test_files/wout_li383_low_res_reference.nc")
DURATION = 2e-3
THREADS = 16


def run_boozer(args):
    case = args.boozer / args.name
    case.mkdir(parents=True, exist_ok=True)
    wout = case / "wout.nc"
    if not wout.exists():
        wout.symlink_to(Path("../wout_li383_low_res_reference.nc"))
    (case / "simple.in").write_text(
        "&config\nnetcdffile='wout.nc'\nisw_field_type=2\n"
        f"integmode={args.integmode}, npoiper2={args.npoiper2}, relerr=1d-12\n"
        f"ntestpart={args.markers}, ntimstep=65, trace_time={args.duration:.17e}\n"
        "facE_al=1000\nstartmode=5, num_surf=2, sbeg=0.1d0, 0.95d0\n"
        "contr_pp=-1d10\ndeterministic=.true.\n/\n")
    env = dict(os.environ, OMP_NUM_THREADS=str(THREADS), OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1", OMP_PROC_BIND="close", OMP_PLACES="cores")
    started = time.monotonic()
    with (case / "simple.log").open("w") as log:
        subprocess.run([str(args.root / "build/simple.x")], cwd=case, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    wall = time.monotonic() - started
    (case / "wall.json").write_text(json.dumps(dict(
        wall_s=wall, integmode=args.integmode, npoiper2=args.npoiper2,
        markers=args.markers, duration_s=args.duration,
        load_average=os.getloadavg())) + "\n")
    print(args.name, wall)


def boozer_case(boozer, name, setup_s):
    meta = json.loads((boozer / name / "wall.json").read_text())
    lost = np.loadtxt(boozer / name / "times_lost.dat")
    t = np.clip(lost[:, 1], 0.0, meta["duration_s"])
    trace = meta["wall_s"] - setup_s
    return dict(model="boozer_gc", integmode=meta["integmode"],
                npoiper2=meta["npoiper2"], markers=meta["markers"],
                losses=int((lost[:, 1] < meta["duration_s"]*(1 - 1e-9)).sum()),
                trace_s=trace, setup_s=setup_s, particle_s=float(t.sum()),
                # Upper bound: assumes all threads busy for the whole trace.
                core_s_per_particle_s=trace*THREADS/float(t.sum()))


def spectre_case(case):
    summary = json.loads((case / "study_summary.json").read_text())
    rows = list(csv.DictReader((case / "study_markers.csv").open(),
                               skipinitialspace=True))
    col = {k: np.array([float(r[k]) for r in rows])
           for k in ("time_s", "wall_s", "fo_fraction", "microsteps")}
    t, w, f, n = col["time_s"], col["wall_s"], col["fo_fraction"], col["microsteps"]
    run = json.loads((case / "run.json").read_text())
    out = dict(model="spectre_" + run["mode"], entry_larmor=run["entry"],
               markers=summary["markers"], losses=summary["physical_losses"],
               numerical_stops=summary["numerical_stops"],
               fo_failures=summary["fo_failures"],
               warning_holds=summary["warning_hold_count"],
               trace_s=summary["trace_s"], setup_s=summary["setup_s"],
               particle_s=float(t.sum()), fo_time_fraction=float((f*t).sum()/t.sum()),
               core_s_per_particle_s=float(w.sum()/t.sum()),
               thread_efficiency=float(w.sum()/(THREADS*summary["trace_s"])),
               load_average=run["load_average"])
    gc = f == 0
    if gc.any():
        out.update(pure_gc_markers=int(gc.sum()),
                   pure_gc_core_s_per_particle_s=float(w[gc].sum()/t[gc].sum()),
                   pure_gc_us_per_step=float(1e6*w[gc].sum()/n[gc].sum()),
                   pure_gc_step_ns=float(1e9*t[gc].sum()/n[gc].sum()))
    return out


def report(args):
    setup = json.loads((args.boozer / "setup_only/wall.json").read_text())["wall_s"]
    cases = {name: boozer_case(args.boozer, name, setup)
             for name in ("x10_default", "x10_mid256", "x10_matched")}
    cases.update({name: spectre_case(args.spectre / name)
                  for name in ("hybrid", "gcprobe", "fo")})
    ref = cases["x10_default"]["core_s_per_particle_s"]
    same_step = cases["x10_matched"]["core_s_per_particle_s"]
    for c in cases.values():
        c["vs_boozer_default"] = c["core_s_per_particle_s"]/ref
        c["vs_boozer_same_step"] = c["core_s_per_particle_s"]/same_step
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(json.dumps(dict(
        wout_url=WOUT_URL, duration_s=DURATION, threads=THREADS,
        cases=cases), indent=1) + "\n")
    for name, c in cases.items():
        print(f"{name:12s} {c['core_s_per_particle_s']:10.1f} "
              f"{c['vs_boozer_default']:7.1f}x {c['vs_boozer_same_step']:6.1f}x")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run")
    r.add_argument("name")
    r.add_argument("integmode", type=int)
    r.add_argument("npoiper2", type=int)
    r.add_argument("markers", type=int)
    r.add_argument("--duration", type=float, default=DURATION)
    r.add_argument("--root", type=Path, default=Path.cwd())
    r.add_argument("--boozer", type=Path, default=Path("build/study/boozer"))
    p = sub.add_parser("report")
    p.add_argument("spectre", type=Path)
    p.add_argument("boozer", type=Path)
    p.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.command == "run":
        args.root = args.root.resolve()
        args.boozer = args.boozer.resolve()
        run_boozer(args)
    else:
        report(args)


if __name__ == "__main__":
    main()
