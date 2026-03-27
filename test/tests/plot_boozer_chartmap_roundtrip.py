#!/usr/bin/env python3
"""Plot Boozer chartmap comparison: fields and orbits side by side."""

import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_field_comparison(prefix, direct_label, chartmap_label, title):
    data = np.loadtxt(f"{prefix}_field_comparison.dat")
    s = data[:, 0]
    theta = data[:, 1]
    Bmod_ref = data[:, 3]
    Bmod_new = data[:, 4]
    rel_err = data[:, 5]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    ax = axes[0]
    ax.scatter(s, Bmod_ref, s=15, label=direct_label, alpha=0.7)
    ax.scatter(s, Bmod_new, s=15, marker="x", label=chartmap_label, alpha=0.7)
    ax.set_xlabel("s")
    ax.set_ylabel("|B| [G]")
    ax.set_title("Bmod vs s")
    ax.legend()

    ax = axes[1]
    ax.scatter(theta, Bmod_ref, s=15, label=direct_label, alpha=0.7)
    ax.scatter(theta, Bmod_new, s=15, marker="x", label=chartmap_label, alpha=0.7)
    ax.set_xlabel("theta_B")
    ax.set_ylabel("|B| [G]")
    ax.set_title("Bmod vs theta")
    ax.legend()

    ax = axes[2]
    ax.semilogy(s, rel_err, "ko", ms=4)
    ax.set_xlabel("s")
    ax.set_ylabel("Relative error")
    ax.set_title("Bmod roundtrip error")
    ax.axhline(1e-5, color="r", ls="--", label="tolerance")
    ax.legend()

    fig.suptitle(f"{title}: field comparison", fontsize=14)
    fig.tight_layout()
    output = f"{prefix}_field_comparison.png"
    fig.savefig(output, dpi=150)
    print(f"Saved {output}")


def plot_orbit_comparison(prefix, direct_label, chartmap_label, title):
    try:
        direct = np.loadtxt(f"{prefix}_orbit_direct.dat")
        chartmap = np.loadtxt(f"{prefix}_orbit_chartmap.dat")
    except Exception as e:
        print(f"Could not load orbit data: {e}")
        return

    n = min(len(direct), len(chartmap))
    direct = direct[:n]
    chartmap = chartmap[:n]

    t_d, s_d, th_d, ph_d = direct[:, 0], direct[:, 1], direct[:, 2], direct[:, 3]
    t_c, s_c, th_c, ph_c = chartmap[:, 0], chartmap[:, 1], chartmap[:, 2], chartmap[:, 3]
    ds = np.abs(s_d - s_c)
    dtheta = th_c - th_d
    dphi = ph_c - ph_d

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(t_d, s_d, label=direct_label, lw=1)
    ax.plot(t_c, s_c, "--", label=chartmap_label, lw=1)
    ax.set_xlabel("time")
    ax.set_ylabel("s")
    ax.set_title("Flux surface label")
    ax.legend()

    ax = axes[0, 1]
    ax.semilogy(t_d, ds, "k-", lw=0.8)
    ax.axhline(1e-4, color="r", ls="--", label="tolerance")
    ax.set_xlabel("time")
    ax.set_ylabel("|s_direct - s_chartmap|")
    ax.set_title("Flux-surface difference")
    ax.legend()

    ax = axes[1, 0]
    ax.plot(t_d, dtheta, color="C2", lw=0.8)
    ax.axhline(0.0, color="0.5", ls="--", lw=0.8)
    ax.set_xlabel("time")
    ax.set_ylabel("theta_B chartmap - direct")
    ax.set_title("Poloidal gauge offset")

    ax = axes[1, 1]
    ax.plot(t_d, dphi, color="C3", lw=0.8)
    ax.axhline(0.0, color="0.5", ls="--", lw=0.8)
    ax.set_xlabel("time")
    ax.set_ylabel("phi_B chartmap - direct")
    ax.set_title("Toroidal gauge offset")

    fig.suptitle(f"{title}: orbit comparison", fontsize=14)
    fig.tight_layout()
    output = f"{prefix}_orbit_comparison.png"
    fig.savefig(output, dpi=150)
    print(f"Saved {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--prefix",
        default="/tmp/boozer_chartmap_roundtrip",
        help="Artifact prefix used by the Fortran comparison test",
    )
    parser.add_argument("--direct-label", default="Direct (VMEC)")
    parser.add_argument("--chartmap-label", default="Chartmap")
    parser.add_argument("--title", default="Boozer chartmap comparison")
    args = parser.parse_args()

    prefix = str(Path(args.prefix))
    plot_field_comparison(prefix, args.direct_label, args.chartmap_label, args.title)
    plot_orbit_comparison(prefix, args.direct_label, args.chartmap_label, args.title)
