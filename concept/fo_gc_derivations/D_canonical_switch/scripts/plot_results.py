#!/usr/bin/env python3
"""Replot retained local switch accuracy; no production-orbit claim."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

data = Path(__file__).resolve().parents[1] / "results/native"
study = json.loads((data / "switch_study.json").read_text())
pullback = json.loads((data / "pullback.json").read_text())
plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 170})
fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), layout="constrained")
for row, (field, title) in enumerate([("ncsx", "NCSX"), ("g3v3", "G3V3L3Fi")]):
    cases = [c for c in study["cases"] if c["field"].startswith(field)]
    eps = np.asarray(cases[0]["eps"])
    for key, label, color, style in [
            ("order0", "Zeroth-order inverse", "#777777", ":"),
            ("energy_matching", "Energy-matching inverse", "#D55E00", "--"),
            ("order1", "First-order inverse", "#0072B2", "-")]:
        values = np.max([c[key] for c in cases], axis=0)
        axes[row, 0].loglog(eps, values, marker="o", color=color,
                            linestyle=style, label=label)
    axes[row, 0].set_ylabel("Max relative magnetic-moment defect")
    axes[row, 0].set_title(f"{title}: resolved FO over five gyroperiods")
    for order, color, style in [(0, "#777777", ":"), (1, "#0072B2", "-")]:
        records = [r for r in pullback["records"]
                   if r["fixture"] == field and r["order"] == order]
        eps_form = sorted({r["epsilon"] for r in records})
        values = [max(r["defect"] for r in records if r["epsilon"] == e)
                  for e in eps_form]
        axes[row, 1].loglog(eps_form, values, marker="s", color=color,
                            linestyle=style, label=f"Selector {order}")
    axes[row, 1].set_title(f"{title}: independent two-form pullback")
    axes[row, 1].set_ylabel("Max absolute component defect")
    axes[row, 1].set_xticks([0.0025, 0.005, 0.01, 0.02],
                           labels=["0.0025", "0.005", "0.01", "0.02"])
    axes[row, 1].xaxis.set_minor_formatter(NullFormatter())
    for ax in axes[row]:
        ax.set_xlabel("Epsilon in fixed study normalization")
        ax.grid(which="both", alpha=0.2)
        ax.legend(fontsize=8)
fig.suptitle("Local generated-switch checks: maxima over three volume samples\n"
             "First-order defects scale quadratically; production hybrid integration remains open",
             fontsize=12)
fig.savefig(data / "switch_accuracy.png")
