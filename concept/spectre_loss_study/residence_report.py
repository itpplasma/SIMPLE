#!/usr/bin/env python3
"""Summarize hybrid cost, loss agreement and FO residence against an FO reference.

usage: residence_report.py OUTPUT REFERENCE CASE [CASE ...]

Each argument is a run.py output directory with matched starts. Cost is
core-seconds per simulated particle-second; loss agreement is the paired status
difference and the maximum difference of the cumulative loss fraction. Cases
with study_episodes.csv also get the FO episode breakdown.
"""
import csv
import json
from pathlib import Path
import sys

import numpy as np

DURATION_BINS_S = (0.0, 1e-6, 1e-5, 1e-4, 1e-3, np.inf)


def rows(path):
    return list(csv.DictReader(path.open(), skipinitialspace=True))


def markers(case):
    r = rows(case / "study_markers.csv")
    return {k: np.array([float(x[k]) for x in r])
            for k in ("time_s", "status", "fo_fraction", "wall_s", "microsteps")}


def loss_curve(m, duration, grid):
    lost = m["time_s"] < duration*(1 - 1e-9)
    times = m["time_s"][lost]
    return lost, np.array([(times <= g).sum() for g in grid])/len(lost)


def episodes(case, mvol):
    r = rows(case / "study_episodes.csv")
    fo_time = np.array([float(x["fo_time_s"]) for x in r])
    rho = np.array([float(x["rho_enter"]) for x in r])
    kind = np.array([int(float(x["end_kind"])) for x in r])
    near = np.where(rho < 0, -1, np.rint(rho).astype(int))
    total = fo_time.sum()
    names = {-1: "in_step", 0: "axis", mvol: "wall"}
    by_place = {names.get(i, f"interface_{i}"): dict(
        count=int((near == i).sum()), fo_time_share=float(fo_time[near == i].sum()/total),
        median_us=float(np.median(fo_time[near == i])*1e6))
        for i in sorted(set(near))}
    by_end = {name: dict(count=int((kind == k).sum()),
                         fo_time_share=float(fo_time[kind == k].sum()/total))
              for k, name in ((1, "exit"), (2, "loss"), (3, "run_end"), (4, "in_step"))}
    by_duration = [dict(lower_s=lo, upper_s=None if np.isinf(hi) else hi,
                        count=int(((fo_time >= lo) & (fo_time < hi)).sum()),
                        fo_time_share=float(fo_time[(fo_time >= lo) & (fo_time < hi)].sum()/total))
                   for lo, hi in zip(DURATION_BINS_S[:-1], DURATION_BINS_S[1:])]
    return dict(episodes=len(r), by_place=by_place, by_end=by_end, by_duration=by_duration)


def summarize(case, reference, ref_curve, ref_lost, grid):
    run = json.loads((case / "run.json").read_text())
    summary = json.loads((case / "study_summary.json").read_text())
    m = markers(case)
    duration = summary["duration_s"]
    lost, curve = loss_curve(m, duration, grid)
    particle_s = m["time_s"].sum()
    fo_share = m["fo_fraction"]*m["time_s"]
    ordered = np.cumsum(np.sort(fo_share)[::-1])/max(fo_share.sum(), np.finfo(float).tiny)
    out = dict(mode=run["mode"], entry_larmor=run["entry"], npoiper2=run["npoiper2"],
               integmode=run.get("integmode", 3), markers=summary["markers"],
               losses=summary["physical_losses"],
               numerical_stops=summary["numerical_stops"],
               fo_failures=summary["fo_failures"],
               warning_holds=summary["warning_hold_count"],
               fo_entries=summary["fo_entries"], trace_s=summary["trace_s"],
               load_average=run["load_average"],
               fo_time_fraction=float(fo_share.sum()/particle_s),
               core_s_per_particle_s=float(m["wall_s"].sum()/particle_s),
               markers_over_90pct_fo=int((m["fo_fraction"] > 0.9).sum()),
               markers_holding_half_fo_time=int(np.searchsorted(ordered, 0.5)) + 1)
    if case != reference:
        out.update(status_differences=int((lost != ref_lost).sum()),
                   loss_curve_max_difference=float(np.abs(curve - ref_curve).max()))
    gc = m["fo_fraction"] == 0
    if gc.any():
        out["pure_gc_core_s_per_particle_s"] = float(m["wall_s"][gc].sum()/m["time_s"][gc].sum())
    if (case / "study_episodes.csv").exists() and out["mode"] == "hybrid":
        out["fo_episodes"] = episodes(case, mvol=3)
    return out


def main():
    output, reference, *cases = map(Path, sys.argv[1:])
    duration = json.loads((reference / "study_summary.json").read_text())["duration_s"]
    grid = np.linspace(0.0, duration, 401)
    ref_lost, ref_curve = loss_curve(markers(reference), duration, grid)
    result = {c.name: summarize(c, reference, ref_curve, ref_lost, grid)
              for c in [reference, *cases]}
    ref_cost = result[reference.name]["core_s_per_particle_s"]
    for c in result.values():
        c["speedup_vs_reference"] = ref_cost/c["core_s_per_particle_s"]
    output.mkdir(parents=True, exist_ok=True)
    (output / f"summary_{reference.name}.json").write_text(
        json.dumps(dict(reference=reference.name, cases=result), indent=1) + "\n")
    for name, c in result.items():
        print(f"{name:12s} loss {c['losses']:4d} dstat {c.get('status_differences', 0):3d} "
              f"curve {c.get('loss_curve_max_difference', 0):.3f} "
              f"FO {c['fo_time_fraction']:6.1%} cost {c['core_s_per_particle_s']:7.0f} "
              f"x{c['speedup_vs_reference']:5.1f}")


if __name__ == "__main__":
    main()
