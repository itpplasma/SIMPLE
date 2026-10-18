#!/usr/bin/env python3
"""Render documentation figures from their versioned text sources."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from make_front import fetch_wout, plot_results

ROOT = Path(__file__).resolve().parents[4]
FIGURES = ROOT / "docs/_static/figures"


def load_front_data():
    source = FIGURES / "data/front"
    metadata = json.loads((source / "metadata.json").read_text())
    losses = np.atleast_1d(
        np.genfromtxt(source / "losses.csv", delimiter=",", names=True)
    )
    confinement = np.loadtxt(source / "confinement.csv", delimiter=",", skiprows=1)
    count = metadata["markers"]
    indices = losses["marker"].astype(int)
    if (
        metadata["numerical_failures"] != 0
        or len(indices) != metadata["losses"]
        or len(np.unique(indices)) != len(indices)
        or np.any(indices < 0)
        or np.any(indices >= count)
    ):
        raise ValueError("Invalid reference marker population")
    trapped = losses["trapped"].astype(bool)
    remaining = np.setdiff1d(np.arange(count), indices)
    passing_remaining = metadata["passing_markers"] - np.count_nonzero(~trapped)
    if not 0 <= passing_remaining <= len(remaining):
        raise ValueError("Inconsistent passing-marker count")
    # Only class counts are needed for non-lost markers, not their orbits.
    classes = np.ones(count)
    classes[indices] = np.where(trapped, 1, -1)
    classes[remaining[:passing_remaining]] = -1
    exit_codes = np.zeros(count, dtype=int)
    exit_codes[indices] = 1
    times = np.full(count, confinement[-1, 0])
    times[indices] = losses["time_s"]
    zend = np.zeros((count, 5))
    xyz = np.zeros((count, 3))
    for column, name in enumerate(("s", "theta_rad", "phi_rad")):
        zend[indices, column] = losses[name]
    for column, name in enumerate(("x_cm", "y_cm", "z_cm")):
        xyz[indices, column] = losses[name]
    return {
        "zend": zend,
        "xend_cart": xyz,
        "times_lost": times,
        "orbit_exit_code": exit_codes,
        "trap_par": classes,
        "version": metadata["plot_label"],
    }, confinement


def main():
    subprocess.run(
        [sys.executable, str(FIGURES / "scripts/make_pipeline.py")], check=True
    )
    work = ROOT / "build/figures/front"
    work.mkdir(parents=True, exist_ok=True)
    fetch_wout(work / "wout.nc")
    results, confinement = load_front_data()
    plot_results(results, confinement, work / "wout.nc", FIGURES / "front.png")
    performance = ROOT / "docs/performance/efss-hpc"
    if (performance / "strong-scaling.csv").exists():
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools/perf/plot_efss_scaling.py"),
                str(performance / "strong-scaling.csv"),
                str(performance / "strong-scaling.png"),
            ],
            check=True,
        )
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools/perf/plot_efss_throughput.py"),
                "--cpu",
                str(performance / "markers64.csv"),
                str(performance / "markers128.csv"),
                str(performance / "strong-scaling.csv"),
                "--output",
                str(performance / "marker-throughput.png"),
            ],
            check=True,
        )


if __name__ == "__main__":
    main()
