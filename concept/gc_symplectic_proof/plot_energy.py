#!/usr/bin/env python3
"""Replot the retained NCSX energy and domain controls."""
from pathlib import Path
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

folder = Path(__file__).with_name("results") / "energy"
coarse = np.genfromtxt(folder / "ncsx_boozer.csv", delimiter=",", names=True,
                       dtype=None, encoding="utf-8")
fine = np.genfromtxt(folder / "half_ncsx_boozer.csv", delimiter=",", names=True,
                     dtype=None, encoding="utf-8")
base_h = coarse["h"][0]
records = json.loads((folder / "summary.json").read_text())["records"][
    "with_domain_and_step_control"]
plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 170})
fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
for column, (data, title) in enumerate([(coarse, "Original fixed step h"),
                                       (fine, "Half step; same elapsed interval")]):
    passing = data[(data["start"] == 1) & (data["method"] == "GAUSS2")]
    time = passing["step"] * passing["h"] / base_h / 1e3
    axes[0, column].plot(time, passing["relative_energy"], color="#0072B2")
    axes[0, column].set_ylabel("Passing start: relative energy error")
    axes[0, column].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[0, column].set_title(title)
    for method, color in [("MIDPOINT", "#D55E00"), ("GAUSS2", "#0072B2"),
                          ("GAUSS4", "#009E73")]:
        rows = data[(data["start"] == 2) & (data["method"] == method)]
        axes[1, column].plot(rows["step"]*rows["h"]/base_h/1e3, rows["r"],
                             color=color, label=method)
        prefix = "81: " if column == 0 else "84: "
        pattern = (re.escape(prefix) + r"ENERGY " + method +
                   r"\s+refused at step (\d+) status=\d+ last radial seed=\s*(\S+)")
        for record in records:
            match = re.fullmatch(pattern, record)
            if match:
                stop_time = int(match[1])*(1 if column == 0 else 0.5)/1e3
                axes[1, column].scatter(stop_time, float(match[2]), color=color,
                                        marker="x", zorder=4)
    axes[1, column].axhline(1, color="#333333", linestyle=":", label="Chart edge")
    axes[1, column].set_ylabel("Trapped start: radial coordinate s")
    axes[1, column].set_ylim(0, 1.05)
    axes[1, column].legend(fontsize=8)
    for row in range(2):
        axes[row, column].set_xlim(0, 100)
        axes[row, column].set_xlabel("Elapsed interval / base step [1000 steps]")
        axes[row, column].grid(alpha=0.2)
fig.suptitle("NCSX Boozer GC controls: finite step and smooth-domain limits\n"
             "Sparse samples; per-step extrema retained. Crosses: boundary/event stops.",
             fontsize=12)
fig.savefig(folder / "energy_domain.png")
