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
if (args.input / "traces.npz").exists():
    retained = np.load(args.input / "traces.npz", allow_pickle=False)
    for name, case in json.loads((args.input / "summary.json").read_text()).items():
        cases[name] = dict(summary=case, **{key: retained[name + "_" + key]
                                          for key in ["trace", "markers", "starts"]})
else:
    for path in sorted(args.input.glob("*/study_summary.json")):
        case = json.loads(path.read_text())
        metadata = json.loads((path.parent / "run.json").read_text())
        arrays = {key: np.genfromtxt(path.parent / f"study_{key}.csv",
                                    delimiter=",", names=True)
                  for key in ["trace", "markers", "starts"]}
        markers = arrays["markers"]
        first_rows = arrays["trace"][np.r_[True, np.diff(arrays["trace"]["marker"]) != 0]]
        gc_rows = first_rows["fo"] == 0
        initial_error = []
        for coordinate in ["rho", "theta", "zeta"]:
            delta = first_rows[coordinate][gc_rows] - arrays["starts"][coordinate][gc_rows]
            if coordinate != "rho":
                delta = np.angle(np.exp(1j*delta))
            initial_error.append(float(np.max(np.abs(delta))) if len(delta) else None)
        case.update(run=metadata, mean_fo_fraction=float(np.mean(markers["fo_fraction"])),
                    fo_time_fraction=float(np.sum(markers["fo_fraction"]*markers["time_s"])
                                           /np.sum(markers["time_s"])),
                    maximum_shell_error=float(np.max(markers["shell_error"])),
                    initial_gc_reference_error_rho_theta_zeta=initial_error,
                    initial_fo_launches=(case["markers"] if metadata["mode"] == "fo" else 0))
        cases[path.parent.name] = dict(summary=case, **arrays)
for name, data in cases.items():
    variant = name.rsplit("_", 1)[-1]
    raw = cases[name.rsplit("_", 1)[0] + "_raw"]["trace"]
    raw_markers = cases[name.rsplit("_", 1)[0] + "_raw"]["markers"]
    markers = data["markers"]
    trace = data["trace"]
    # Compare only identical marker/time/state-representation records. A mode
    # switch changes particle vs guiding-centre positions, so exclude those
    # pairs rather than treating a gyroradius offset as interpolation error.
    keys = lambda t: zip(t["marker"].astype(int), np.round(t["time_s"], 15))
    lookup = dict(zip(keys(raw), range(len(raw))))
    distances = []
    endpoint = []
    representation_mismatches = 0
    unmatched_records = 0
    for j, key in enumerate(keys(trace)):
        i = lookup.get(key)
        if i is None:
            unmatched_records += 1
            continue
        if trace["fo"][j] != raw["fo"][i]:
            representation_mismatches += 1
            continue
        distance = np.sqrt(sum((trace[k][j] - raw[k][i])**2
                               for k in ["x_cm", "y_cm", "z_cm"]))
        distances.append(distance)
        if j == len(trace) - 1 or trace["marker"][j + 1] != trace["marker"][j]:
            endpoint.append(distance)
    data["summary"].update(maximum_paired_position_cm=max(distances, default=None),
                           median_endpoint_position_cm=(float(np.median(endpoint))
                                                        if endpoint else None),
                           paired_records=len(distances),
                           representation_mismatches=representation_mismatches,
                           unmatched_records=unmatched_records)
    both_lost = (markers["status"] == 88) & (raw_markers["status"] == 88)
    data["summary"].update(
        paired_status_disagreements=int(np.count_nonzero(
            markers["status"] != raw_markers["status"])),
        maximum_paired_loss_time_s=(float(np.max(np.abs(
            markers["time_s"][both_lost] - raw_markers["time_s"][both_lost])))
            if both_lost.any() else None))
summaries = {name: data["summary"] for name, data in cases.items()}
(args.output / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
np.savez_compressed(args.output / "traces.npz", **{
    name + "_" + key: data[key] for name, data in cases.items()
    for key in ["trace", "markers", "starts"]})
styles = {"raw": ("#595959", "main / spectral", "", None),
          "coarse": ("#D55E00", "16×32×16", "//", 16),
          "spline": ("#0072B2", "32×64×32", "..", 32),
          "fine": ("#009E73", "48×96×48", "xx", 48)}
variants = [v for v in styles if f"ncsx_sheet_fo_{v}" in summaries]
plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 170})
fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
for row, field in enumerate(["ncsx", "g3v3"]):
    for col, sampling in enumerate(["sheet", "volume"]):
        ax = axes[row, col]
        for index, variant in enumerate(variants):
            values = [summaries[f"{field}_{sampling}_{mode}_{variant}"]["trace_s"]
                      for mode in ["fo", "hybrid"]]
            ax.bar(np.arange(2) + (index - (len(variants) - 1)/2)*0.18,
                   values, 0.17, color=styles[variant][0],
                   hatch=styles[variant][2], label=styles[variant][1])
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
    grid_sizes = [styles[v][3] for v in variants[1:]]
    axes[0].semilogy(grid_sizes,
                     [v["maximum_paired_position_cm"]*1e4 for v in values],
                     marker=["o", "s", "^", "D"][index], label=label)
    axes[1].semilogy(grid_sizes, [v["maximum_shell_error"] for v in values],
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
fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
for row, field in enumerate(["ncsx", "g3v3"]):
    for col, sampling in enumerate(["sheet", "volume"]):
        ax = axes[row, col]
        for variant in [v for v in variants if v != "coarse"]:
            for mode in ["fo", "hybrid"]:
                name = f"{field}_{sampling}_{mode}_{variant}"
                markers = cases[name]["markers"]
                lost = np.sort(markers["time_s"][markers["status"] == 88])
                times = np.r_[0, lost, summaries[name]["duration_s"]]*1e6
                fraction = np.r_[0, np.arange(1, len(lost) + 1), len(lost)]/len(markers)
                ax.step(times, fraction, where="post", color=styles[variant][0],
                        linestyle="-" if mode == "fo" else "--",
                        label=f"{styles[variant][1]} {'FO' if mode == 'fo' else 'hybrid'}"
                        f" ({len(lost)}/{len(markers)})")
        ax.set_title(f"{'NCSX' if field == 'ncsx' else 'G3V3L3Fi'} · {sampling} starts")
        ax.set_xlabel("Physical time [µs]")
        ax.set_ylabel("Cumulative physical loss / initial markers")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7.5)
fig.suptitle("Physical outer-edge losses · matched GC starts within each sampling set\n"
             "Solid: explicit FO reference; dashed: current production hybrid", fontsize=13)
fig.savefig(args.output / "losses.png")
plt.close(fig)
for name, case in summaries.items():
    print(name, f"trace={case['trace_s']:.4f}s",
          f"FO={case['mean_fo_fraction']:.3f}",
          f"loss/stop/fail/hold={case['physical_losses']}/{case['numerical_stops']}/"
          f"{case['fo_failures']}/{case['warning_hold_count']}")
