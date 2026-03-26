#!/usr/bin/env python3
"""Plot Boozer chartmap roundtrip comparison: fields and orbits side by side."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_field_comparison():
    data = np.loadtxt("/tmp/boozer_field_comparison.dat")
    s = data[:, 0]
    theta = data[:, 1]
    Bmod_ref = data[:, 3]
    Bmod_new = data[:, 4]
    rel_err = data[:, 5]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    ax = axes[0]
    ax.scatter(s, Bmod_ref, s=15, label="Direct (VMEC)", alpha=0.7)
    ax.scatter(s, Bmod_new, s=15, marker="x", label="Chartmap", alpha=0.7)
    ax.set_xlabel("s")
    ax.set_ylabel("|B| [G]")
    ax.set_title("Bmod vs s")
    ax.legend()

    ax = axes[1]
    ax.scatter(theta, Bmod_ref, s=15, label="Direct (VMEC)", alpha=0.7)
    ax.scatter(theta, Bmod_new, s=15, marker="x", label="Chartmap", alpha=0.7)
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

    fig.suptitle("Boozer Chartmap Field Roundtrip", fontsize=14)
    fig.tight_layout()
    fig.savefig("/tmp/boozer_field_comparison.png", dpi=150)
    print("Saved /tmp/boozer_field_comparison.png")


def plot_orbit_comparison():
    try:
        direct = np.loadtxt("/tmp/orbit_direct.dat")
        chartmap = np.loadtxt("/tmp/orbit_chartmap.dat")
    except Exception as e:
        print(f"Could not load orbit data: {e}")
        return

    n = min(len(direct), len(chartmap))
    direct = direct[:n]
    chartmap = chartmap[:n]

    t_d, s_d, th_d, ph_d = direct[:, 0], direct[:, 1], direct[:, 2], direct[:, 3]
    t_c, s_c, th_c, ph_c = chartmap[:, 0], chartmap[:, 1], chartmap[:, 2], chartmap[:, 3]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(t_d, s_d, label="Direct", lw=1)
    ax.plot(t_c, s_c, "--", label="Chartmap", lw=1)
    ax.set_xlabel("time")
    ax.set_ylabel("s")
    ax.set_title("Flux surface label")
    ax.legend()

    ax = axes[0, 1]
    ax.plot(t_d, th_d, label="Direct", lw=0.5)
    ax.plot(t_c, th_c, "--", label="Chartmap", lw=0.5)
    ax.set_xlabel("time")
    ax.set_ylabel("theta_B")
    ax.set_title("Poloidal angle")
    ax.legend()

    ax = axes[1, 0]
    ax.plot(t_d, ph_d, label="Direct", lw=0.5)
    ax.plot(t_c, ph_c, "--", label="Chartmap", lw=0.5)
    ax.set_xlabel("time")
    ax.set_ylabel("phi_B")
    ax.set_title("Toroidal angle")
    ax.legend()

    ax = axes[1, 1]
    ax.semilogy(t_d, np.abs(s_d - s_c), "k-", lw=0.5)
    ax.set_xlabel("time")
    ax.set_ylabel("|s_direct - s_chartmap|")
    ax.set_title("Orbit difference")
    ax.axhline(1e-4, color="r", ls="--", label="tolerance")
    ax.legend()

    fig.suptitle("Boozer Chartmap Orbit Roundtrip", fontsize=14)
    fig.tight_layout()
    fig.savefig("/tmp/boozer_orbit_comparison.png", dpi=150)
    print("Saved /tmp/boozer_orbit_comparison.png")


if __name__ == "__main__":
    plot_field_comparison()
    plot_orbit_comparison()
