#!/usr/bin/env python3
"""Retain physical full-duration controls and plot losses and measured cost."""
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
parser.add_argument("--particle-boundary", action="store_true",
                    help="Hybrid cases include the outer-boundary FO repair")
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
cases, arrays = {}, {}
retained = None
if (args.input / "traces.npz").exists():
    retained = json.loads((args.input / "summary.json").read_text())
    with np.load(args.input / "traces.npz", allow_pickle=False) as data:
        arrays = {key: data[key] for key in data.files}
for field in ["ncsx", "g3v3"]:
    for model in ["fo_default", "fo_fine", "hybrid_default"]:
        name = field + "_" + model
        folder = args.input / name
        if retained is None:
            case = json.loads((folder / "study_summary.json").read_text())
            case["run"] = json.loads((folder / "run.json").read_text())
            for key in ["markers", "starts", "trace"]:
                arrays[name + "_" + key] = np.genfromtxt(
                    folder / f"study_{key}.csv", delimiter=",", names=True)
        else:
            case = retained["cases"][name]
        markers = arrays[name + "_markers"]
        case.update(fo_time_fraction=float(np.sum(markers["fo_fraction"]*markers["time_s"])
                                           /np.sum(markers["time_s"])),
                    maximum_shell_error=float(np.max(markers["shell_error"])),
                    worker_wall_fraction=float(np.sum(markers["wall_s"])
                                               /(case["threads"]*case["trace_s"])))
        if any(case[k] for k in ["numerical_stops", "fo_failures", "warning_hold_count"]):
            raise RuntimeError(f"Unresolved physical control: {name}")
        cases[name] = case
particle_boundary = args.particle_boundary or bool(
    retained and retained.get("particle_boundary", False))
comparisons = {}
for field in ["ncsx", "g3v3"]:
    reference = arrays[field + "_fo_fine_markers"]
    for model in ["fo_default", "hybrid_default"]:
        name = field + "_" + model
        markers = arrays[name + "_markers"]
        difference = (markers["status"] == 88).astype(int) - (reference["status"] == 88)
        both = (markers["status"] == 88) & (reference["status"] == 88)
        comparisons[name] = dict(
            reference=field + "_fo_fine", paired_mismatches=int(np.count_nonzero(difference)),
            paired_loss_fraction_difference=float(np.mean(difference)),
            paired_standard_error=float(np.std(difference, ddof=1)/np.sqrt(len(markers))),
            maximum_paired_loss_time_s=(float(np.max(np.abs(
                markers["time_s"][both] - reference["time_s"][both]))) if both.any() else None),
            zero_mismatch_95_percent_upper_bound=(1 - 0.05**(1/len(markers))
                                                  if not difference.any() else None))
boundary_scope = ("FO entry resolves the particle absorbing outer edge in the hybrid."
                  if particle_boundary else "FO stops at the particle edge; GC "
                  "segments use the guiding-centre edge convention.")
summary = dict(cases=cases, comparisons=comparisons, particle_boundary=particle_boundary,
               scope="Current hybrid controls, not a phase-consistent or symplectic proof; "
                     "volume-sampled GC starts with 0.1 <= rho <= Mvol-0.05 and the "
                     "deterministic production entry phase. " + boundary_scope)
(args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
np.savez_compressed(args.output / "traces.npz", **arrays)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 170})
styles = [("fo_default", "FO, default grid", "#0072B2", ":"),
          ("fo_fine", "FO, refined grid", "#333333", "-"),
          ("hybrid_default", "Current hybrid", "#D55E00", "--")]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
for ax, field in zip(axes, ["ncsx", "g3v3"]):
    for model, label, color, line in styles:
        name = field + "_" + model
        markers = arrays[name + "_markers"]
        lost = np.sort(markers["time_s"][markers["status"] == 88])
        time = np.r_[0, lost, cases[name]["duration_s"]]*1e3
        fraction = np.r_[0, np.arange(1, len(lost) + 1), len(lost)]/len(markers)
        ax.step(time, fraction, where="post", color=color, linestyle=line,
                label=f"{label}: {len(lost)}/{len(markers)}", linewidth=1.8)
        if model == "fo_fine":
            z = 1.959963984540054
            n = len(markers)
            center = (fraction + z*z/(2*n))/(1 + z*z/n)
            half = z*np.sqrt(fraction*(1 - fraction)/n + z*z/(4*n*n))/(1 + z*z/n)
            ax.fill_between(time, center - half, center + half, step="post",
                            color=color, alpha=0.12, label="Pointwise 95% Wilson band")
    ax.set_title("NCSX" if field == "ncsx" else "G3V3L3Fi")
    ax.set_xlabel("Physical time [ms]")
    ax.set_ylabel("Reported edge loss / initial markers")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
first = next(iter(cases.values()))
fig.suptitle(f"3.5 keV · {first['markers']} matched volume-sampled GC starts\n"
             + ("Current switch and entry phase; particle absorbing outer edge"
                if particle_boundary else
                "Current switch and entry phase; GC-centre / FO-particle edge conventions"),
             fontsize=12)
fig.savefig(args.output / "losses.png")
plt.close(fig)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
for ax, field in zip(axes, ["ncsx", "g3v3"]):
    for j, (model, label, _, _) in enumerate(styles):
        case = cases[field + "_" + model]
        overhead = max(0, case["run"]["process_s"] - case["setup_s"] - case["trace_s"])
        left = 0
        for width, color, hatch, title in [
                (overhead, "#999999", "//", "Native launch / sampling / output"),
                (case["setup_s"], "#E69F00", "..", "Field initialization"),
                (case["trace_s"], "#0072B2", "", "Parallel tracing")]:
            ax.barh(j, width, left=left, color=color, hatch=hatch,
                    label=title if j == 0 else None)
            left += width
        ax.text(left, j, f" {left:.1f} s", va="center", fontsize=9)
    ax.set_yticks(range(3), [s[1] for s in styles])
    ax.invert_yaxis()
    ax.set_xlabel("Measured end-to-end wall time [s]")
    ax.set_xlim(0, 1.18*max(case["run"]["process_s"] for case in cases.values()))
    ax.set_title("NCSX" if field == "ncsx" else "G3V3L3Fi")
    ax.grid(axis="x", alpha=0.2)
axes[0].legend(fontsize=8, loc="lower right")
fig.suptitle(f"{first['markers']} markers · {first['duration_s']*1e3:g} ms · "
             f"{first['threads']} OpenMP threads on local shared host", fontsize=13)
fig.savefig(args.output / "cost.png")
plt.close(fig)
for name, case in cases.items():
    print(name, f"trace={case['trace_s']:.2f}s total={case['run']['process_s']:.2f}s",
          f"loss={case['physical_losses']}/{case['markers']} FO={case['fo_time_fraction']:.3f}")
