#!/usr/bin/env python3
"""Retain the five-marker loss diagnosis and render its measured controls."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def compact_trace(trace):
    """Thin plotting records; retain extrema, endpoints and both switch sides."""
    keep = np.zeros(len(trace), dtype=bool)
    for marker in np.unique(trace["marker"]):
        ix = np.flatnonzero(trace["marker"] == marker)
        keep[ix[::16]] = True
        keep[ix[[0, -1]]] = True
        switches = np.flatnonzero(np.diff(trace["fo"][ix]))
        keep[ix[switches]] = True
        keep[ix[switches + 1]] = True
        radial_change = np.diff(trace["rho"][ix])
        extrema = np.flatnonzero(radial_change[1:]*radial_change[:-1] <= 0) + 1
        keep[ix[extrema]] = True
    return trace[keep]


def collect(args):
    previous = args.previous
    with np.load(previous / "traces.npz", allow_pickle=False) as old:
        arrays = {"before_markers": old["g3v3_hybrid_default_markers"],
                  "fo_markers": old["g3v3_fo_fine_markers"]}
    cases = {}
    names = ["fo_phases129", "hybrid_five_extended", "hybrid129_halfstep",
             "hybrid129_quarterstep", "hybrid129_eighthstep", "hybrid129_fine",
             "hybrid129_construction", "hybrid129_wide", "hybrid129_allfo"]
    for name in names:
        folder = args.input / name
        case = json.loads((folder / "study_summary.json").read_text())
        case["run"] = json.loads((folder / "run.json").read_text())
        if any(case[key] for key in
               ["numerical_stops", "fo_failures", "warning_hold_count"]):
            raise RuntimeError(f"Unresolved control: {name}")
        cases[name] = case
        for kind in ["starts", "markers", "trace"]:
            data = np.atleast_1d(np.genfromtxt(
                folder / f"study_{kind}.csv", delimiter=",", names=True))
            arrays[name + "_" + kind] = (compact_trace(data)
                                          if kind == "trace" else data)
    for label in ["default", "fine"]:
        name = "g3v3_hybrid_" + label
        folder = args.wall / name
        case = json.loads((folder / "study_summary.json").read_text())
        case["run"] = json.loads((folder / "run.json").read_text())
        if any(case[key] for key in
               ["numerical_stops", "fo_failures", "warning_hold_count"]):
            raise RuntimeError(f"Unresolved control: {name}")
        cases[name] = case
        arrays[name + "_markers"] = np.genfromtxt(
            folder / "study_markers.csv", delimiter=",", names=True)
    before = arrays["before_markers"]
    fo = arrays["fo_markers"]
    wall = arrays["g3v3_hybrid_default_markers"]
    ids = fo["marker"][before["status"] != fo["status"]].astype(int)
    markers = []
    for marker in ids:
        i = marker - 1
        markers.append(dict(marker=int(marker), fo_loss_s=float(fo[i]["time_s"]),
                            before_status=int(before[i]["status"]),
                            wall_status=int(wall[i]["status"]),
                            wall_time_s=float(wall[i]["time_s"])))
    phases = arrays["fo_phases129_markers"]
    summary = dict(
        previous=str(previous), cases=cases, differing_markers=markers,
        subset_original_ids=[62, 124, 126, 129, 247],
        phase_marker_original_id=129,
        phase_offsets_rad=(2*np.pi*np.arange(len(phases))/len(phases)).tolist(),
        phase_losses_by_2ms=int(np.count_nonzero(phases["time_s"] <= .002)),
        phase_loss_time_range_s=[float(min(phases["time_s"])),
                                 float(max(phases["time_s"]))],
        allfo_minus_phase_zero_loss_s=float(
            arrays["hybrid129_allfo_markers"][0]["time_s"] - phases[0]["time_s"]),
        wall_vs_fo_disagreements=int(np.count_nonzero(wall["status"] != fo["status"])),
        wall_default_vs_fine_disagreements=int(np.count_nonzero(
            wall["status"] != arrays["g3v3_hybrid_fine_markers"]["status"])),
        trace_retention="Every 16th record plus radial extrema, endpoints and GC/FO switches; "
                        "all marker outcomes and starts retained. Raw dense CSVs "
                        "remain under ignored build/study/spectre-loss-diagnosis.",
        conclusion="Four missed particle-wall events are a confirmed boundary "
                   "inconsistency. Marker 129 remains sensitive to the present "
                   "switch/phase model. Its path is not established converged "
                   "by the step controls; do not call the residual a validated "
                   "physical loss difference or a whole-hybrid acceptance.")
    return summary, arrays


def plot(summary, arrays, output):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 180})
    fig = plt.figure(figsize=(11, 8), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[.9, 1])
    ax = fig.add_subplot(grid[0, :])
    rows = summary["differing_markers"]
    y = np.arange(len(rows))
    ax.scatter([r["fo_loss_s"]*1e3 for r in rows], y + .14, color="#222222",
               marker="x", s=45, label="FO reference: particle-wall loss")
    ax.scatter(np.full(len(rows), 2), y - .14, edgecolors="#777777",
               facecolors="none", marker=">", s=45,
               label="Before repair: survives to cutoff")
    for status, fill, label in [(88, "#D55E00", "Wall repair: particle-wall loss"),
                                 (0, "none", "Wall repair: survives to cutoff")]:
        selected = [(i, r) for i, r in enumerate(rows) if r["wall_status"] == status]
        ax.scatter([r["wall_time_s"]*1e3 for _, r in selected],
                   [i for i, _ in selected], marker="s", s=35,
                   edgecolors="#D55E00", facecolors=fill, label=label)
    ax.axvline(2, color="#777777", linestyle=":", linewidth=1)
    ax.set(yticks=y, yticklabels=[str(r["marker"]) for r in rows],
           ylabel="Original marker ID", xlabel="Event or censoring time [ms]",
           xlim=(.35, 2.12), ylim=(-.5, len(rows) - .5),
           title="G3V3L3Fi: four of five missed wall losses recovered (48 → 52 / 256)")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=.2)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.19), ncol=2)

    ax = fig.add_subplot(grid[1, 0])
    phases = arrays["fo_phases129_markers"]
    degrees = np.rad2deg(summary["phase_offsets_rad"])
    for selected, color, marker, label in [
            (phases["time_s"] <= .002, "#0072B2", "o", "Lost by 2 ms"),
            (phases["time_s"] > .002, "#D55E00", "s", "Survives at 2 ms")]:
        ax.scatter(degrees[selected], phases["time_s"][selected]*1e3,
                   color=color, marker=marker, label=label)
    ax.axhline(2, color="#777777", linestyle="--", linewidth=1,
               label="2 ms cutoff")
    ax.set(xlabel="FO launch-phase offset [degrees]", ylabel="Particle-wall loss time [ms]",
           title="Marker 129: physical launch-phase sensitivity", xticks=[0, 90, 180, 270, 360])
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(grid[1, 1])
    ref = arrays["fo_phases129_trace"]
    ref = ref[ref["marker"] == 1]
    ax.plot(ref["time_s"]*1e3, ref["rho"], color="#222222", label="FO, phase zero")
    for name, marker, color, label in [
            ("hybrid_five_extended", 4, "#D55E00", "Hybrid, original microstep"),
            ("hybrid129_eighthstep", 1, "#0072B2", "Hybrid, microstep ≈ h/8")]:
        trace = arrays[name + "_trace"]
        trace = trace[trace["marker"] == marker]
        ax.plot(trace["time_s"]*1e3, np.where(trace["fo"] == 1, trace["rho"], np.nan),
                color=color, label=label)
        ax.plot(trace["time_s"]*1e3, np.where(trace["fo"] == 0, trace["rho"], np.nan),
                color=color, linestyle=":")
    ax.axhline(3, color="#777777", linestyle="--", linewidth=1)
    ax.axvline(2, color="#777777", linestyle=":", linewidth=1)
    ax.set(xlim=(1.75, 2.125), ylim=(2.65, 3.015), xlabel="Physical time [ms]",
           ylabel="Reported stacked radius", title="Marker 129: late orbit window")
    ax.text(.02, .03, "Solid: particle · dotted: GC centre", transform=ax.transAxes,
            fontsize=8, bbox=dict(facecolor="white", edgecolor="none", alpha=.8))
    ax.grid(alpha=.2)
    ax.legend(fontsize=8, loc="upper left")
    fig.suptitle("3.5 keV alpha particles · deterministic current switch\n"
                 "Boundary defect repaired; residual switch-model error remains open", fontsize=13)
    fig.savefig(output / "diagnosis.png")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--wall", type=Path, default=Path("build/study/spectre-full-wall-fixed"))
    parser.add_argument("--previous", type=Path,
                        default=Path(__file__).parent / "results/full-gc-fixed")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if (args.input / "diagnosis.json").exists():
        summary = json.loads((args.input / "diagnosis.json").read_text())
        with np.load(args.input / "traces.npz", allow_pickle=False) as data:
            arrays = {key: data[key] for key in data.files}
    else:
        summary, arrays = collect(args)
    (args.output / "diagnosis.json").write_text(json.dumps(summary, indent=2) + "\n")
    np.savez_compressed(args.output / "traces.npz", **arrays)
    plot(summary, arrays, args.output)
    print(json.dumps({key: value for key, value in summary.items() if key != "cases"}, indent=2))


if __name__ == "__main__":
    main()
