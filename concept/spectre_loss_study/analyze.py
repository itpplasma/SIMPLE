#!/usr/bin/env python3
"""Retain paired pilot summaries and render scientific comparison plots."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(__doc__)
parser.add_argument("input", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
cases = {}
for path in sorted(args.input.glob("*/study_summary.json")):
    case = json.loads(path.read_text())
    metadata = json.loads((path.parent / "run.json").read_text())
    markers = np.genfromtxt(path.parent / "study_markers.csv", delimiter=",", names=True)
    trace = np.genfromtxt(path.parent / "study_trace.csv", delimiter=",", names=True)
    case.update(run=metadata, mean_fo_fraction=float(np.mean(markers["fo_fraction"])),
                maximum_shell_error=float(np.max(markers["shell_error"])))
    cases[path.parent.name] = dict(summary=case, trace=trace)
for name, data in cases.items():
    variant = name.rsplit("_", 1)[-1]
    raw = cases[name.rsplit("_", 1)[0] + "_raw"]["trace"]
    trace = data["trace"]
    # Compare only identical marker/time/state-representation records. A mode
    # switch changes particle vs guiding-centre positions, so exclude those
    # pairs rather than treating a gyroradius offset as interpolation error.
    keys = lambda t: zip(t["marker"].astype(int), np.round(t["time_s"], 15))
    lookup = dict(zip(keys(raw), range(len(raw))))
    distances = []
    endpoint = []
    for j, key in enumerate(keys(trace)):
        i = lookup.get(key)
        if i is None or trace["fo"][j] != raw["fo"][i]:
            continue
        distance = np.sqrt(sum((trace[k][j] - raw[k][i])**2
                               for k in ["x_cm", "y_cm", "z_cm"]))
        distances.append(distance)
        if j == len(trace) - 1 or trace["marker"][j + 1] != trace["marker"][j]:
            endpoint.append(distance)
    data["summary"].update(maximum_paired_position_cm=max(distances, default=None),
                           median_endpoint_position_cm=(float(np.median(endpoint))
                                                        if endpoint else None),
                           paired_records=len(distances))
summaries = {name: data["summary"] for name, data in cases.items()}
(args.output / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
colors = ["#595959", "#D55E00", "#0072B2", "#009E73"]
variants = ["raw", "coarse", "spline", "fine"]
labels = ["main / spectral", "16×32×16", "32×64×32", "48×96×48"]
plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 170})
fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
for row, field in enumerate(["ncsx", "g3v3"]):
    for col, sampling in enumerate(["sheet", "volume"]):
        ax = axes[row, col]
        for index, variant in enumerate(variants):
            values = [summaries[f"{field}_{sampling}_{mode}_{variant}"]["trace_s"]
                      for mode in ["fo", "hybrid"]]
            ax.bar(np.arange(2) + (index - 1.5)*0.18, values, 0.17,
                   color=colors[index], label=labels[index])
        ax.set_xticks([0, 1], ["FO throughout", "GC + sheet FO"])
        ax.set_ylabel("Parallel tracing wall time [s]")
        ax.set_title(f"{'NCSX' if field == 'ncsx' else 'G3V3L3Fi'} · {sampling} starts")
        ax.grid(axis="y", alpha=0.2)
axes[0, 0].legend(fontsize=9)
first = next(iter(summaries.values()))
fig.suptitle(f"3.5 keV pilot · {first['markers']} matched markers · "
             f"{first['duration_s']*1e6:g} µs · {first['threads']} threads\n"
             "Tracing only; initialization and sampling excluded", fontsize=13)
fig.savefig(args.output / "throughput.png")
plt.close(fig)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
for index, (field, sampling) in enumerate((f, s) for f in ["ncsx", "g3v3"]
                                        for s in ["sheet", "volume"]):
    values = [summaries[f"{field}_{sampling}_fo_{v}"] for v in variants[1:]]
    label = f"{'NCSX' if field == 'ncsx' else 'G3V3L3Fi'}, {sampling}"
    axes[0].semilogy([16, 32, 48],
                     [v["maximum_paired_position_cm"]*1e4 for v in values],
                     marker=["o", "s", "^", "D"][index], label=label)
    axes[1].semilogy([16, 32, 48], [v["maximum_shell_error"] for v in values],
                     marker=["o", "s", "^", "D"][index], label=label)
axes[0].set_ylabel("Max paired position difference from main [µm]")
axes[1].set_ylabel("Max normalized kinetic-energy error")
for ax in axes:
    ax.set_xlabel("Radial grid points (angular: 2n × n)")
    ax.set_xticks([16, 32, 48])
    ax.grid(alpha=0.2)
axes[0].legend(fontsize=9)
fig.suptitle("FO throughout · independent spectral main reference\n"
             "Fixed RK45 tolerance 10⁻⁹ and 0.1 rad gyro cap", fontsize=13)
fig.savefig(args.output / "convergence.png")
plt.close(fig)
for name, case in summaries.items():
    print(name, f"trace={case['trace_s']:.4f}s",
          f"FO={case['mean_fo_fraction']:.3f}",
          f"loss/stop/fail/hold={case['physical_losses']}/{case['numerical_stops']}/"
          f"{case['fo_failures']}/{case['warning_hold_count']}")
