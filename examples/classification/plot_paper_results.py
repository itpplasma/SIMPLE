#!/usr/bin/env python3
"""
Reproduce classification plots from the JPP 2020 paper:
  Albert, Kasilov, Kernbichler, "Accelerated methods for direct computation
  of fusion alpha particle losses within stellarator optimization",
  J. Plasma Phys. 86, 815860201 (2020).

Uses a QH (quasi-helical) equilibrium from Landreman & Paul 2021 as proxy
for the QH configuration (Drevlak 2018) in the original paper.

Generates:
  - Figure 6 style: loss time vs trapping parameter for s=0.3 and s=0.6
  - Figure 8 style: orbit classification in (theta, v_par/v) at s=0.6
  - Volume classification: (s, J_perp) plot
  - Table 1 style: regular fractions for trapped/passing orbits

Usage:
    python plot_paper_results.py [base_dir]

    base_dir defaults to /tmp/simple_classification
"""

import sys
from pathlib import Path

import numpy as np

# Import pysimple.plotting without requiring compiled backend
_plotting_path = Path(__file__).resolve().parent.parent.parent / "python" / "pysimple" / "plotting.py"
import importlib.util
_spec = importlib.util.spec_from_file_location("pysimple.plotting", _plotting_path)
_plotting = importlib.util.module_from_spec(_spec)
sys.modules["pysimple.plotting"] = _plotting
_spec.loader.exec_module(_plotting)

load_loss_data = _plotting.load_loss_data
plot_kde_loss_density = _plotting.plot_kde_loss_density

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde


def plot_paper_loss_figure(data_s03, data_s06, output_path):
    """
    Figure 5-7 style: two-panel loss plot for s=0.3 and s=0.6.

    Each panel shows:
      - Black dots for lost particles (log10(t) vs trapping parameter)
      - Blue KDE density shading
      - Red confined fraction curve with error bands (right axis)
      - Dashed line at trapped-passing boundary (theta_trap=0)
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, data, label in [(axes[0], data_s03, "s = 0.3"),
                             (axes[1], data_s06, "s = 0.6")]:
        lost = data.lost_mask
        tlost = np.abs(data.loss_times[lost])
        trap_par = data.trap_parameter[lost]
        t_threshold = 0.99 * data.trace_time
        kde_mask = tlost < t_threshold

        # KDE density
        if np.sum(kde_mask) > 10:
            log_tlost = np.log10(tlost[kde_mask])
            trap_kde = trap_par[kde_mask]
            try:
                kde = gaussian_kde([log_tlost, trap_kde])
                X, Y = np.mgrid[-5:0:100j, -2:1:150j]
                positions = np.vstack([X.ravel(), Y.ravel()])
                Z = np.reshape(kde(positions).T, X.shape)
                ax.imshow(
                    np.rot90(Z), cmap=plt.cm.Blues,
                    extent=[-5, 0, -2, 1], aspect="auto")
            except np.linalg.LinAlgError:
                pass

        # Lost particle dots
        if np.sum(lost) > 0:
            ax.plot(np.log10(tlost), trap_par, "k,", markersize=0.5, alpha=0.5)

        ax.axhline(y=0, color="k", linestyle="--", linewidth=0.8)
        ax.set_ylim([-2.0, 1.0])
        ax.set_xlim([-5.0, 0.0])
        ax.set_ylabel(r"trapping parameter $\theta_\mathrm{trap}$")
        ax.set_xlabel(r"$t$ / s")
        ax.set_xticks([-5, -4, -3, -2, -1, 0])
        ax.set_xticklabels(
            [r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$",
             r"$10^{-2}$", r"$10^{-1}$", r"$1$"])
        ax.set_title(label)

        # Confined fraction on twin axis
        ax2 = ax.twinx()
        n = data.n_particles
        t_grid = np.logspace(-5, 0, 200)
        fc = np.array([np.sum(np.abs(data.loss_times) >= t) / n for t in t_grid])
        sigma = np.sqrt(fc * (1 - fc) / n)

        ax2.plot(np.log10(t_grid), fc, color="tab:red", lw=1.5)
        ax2.fill_between(
            np.log10(t_grid), fc - 1.96 * sigma, fc + 1.96 * sigma,
            color="tab:red", alpha=0.15)
        ax2.set_ylim([0.6, 1.05])
        ax2.set_ylabel("confined fraction $f_c$", color="tab:red")
        ax2.tick_params("y", colors="tab:red")

    fig.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_classification_theta_pitch(data_dir, output_path):
    """
    Figure 8 style: orbit classification in (theta, v_par/v) space at s=0.6.

    Reads class_parts.dat and start.dat from a single-surface run.
    Colors:
      - Green background: regular orbits (filled square)
      - Blue circle: early losses (before tcut)
      - Red cross: chaotic orbits with late losses
      - White line: trapped-passing boundary
    """
    data_dir = Path(data_dir)

    class_data = np.loadtxt(data_dir / "class_parts.dat")
    iclass_jpar = class_data[:, 3].astype(int)

    tl_data = np.loadtxt(data_dir / "times_lost.dat")
    loss_times = tl_data[:, 1]
    trap_par = tl_data[:, 2]

    start_file = data_dir / "start.dat"
    if start_file.exists():
        start_data = np.loadtxt(start_file)
        theta = np.mod(start_data[:, 1], 2 * np.pi) / np.pi
        pitch = start_data[:, 4]
    else:
        print("  WARNING: start.dat not found, skipping Figure 8 plot")
        return

    # Read trace_time from simple.in
    trace_time = 1.0
    simple_in = data_dir / "simple.in"
    if simple_in.exists():
        for line in simple_in.read_text().splitlines():
            if "trace_time" in line and "=" in line:
                val = line.split("=")[1].split("!")[0].strip().rstrip(",")
                trace_time = float(val.replace("d", "e"))

    # Read tcut from simple.in
    tcut = 0.1
    if simple_in.exists():
        for line in simple_in.read_text().splitlines():
            if "tcut" in line and "=" in line and "ntcut" not in line:
                val = line.split("=")[1].split("!")[0].strip().rstrip(",")
                tcut = float(val.replace("d", "e"))

    fig, ax = plt.subplots(figsize=(6, 5))

    # Classify particles into categories matching the paper
    regular = iclass_jpar == 1
    prompt_loss = iclass_jpar == 0
    chaotic = iclass_jpar == 2

    # For chaotic particles, split into early and late losses
    lost_early = (loss_times > 0) & (loss_times < tcut)
    lost_late = (loss_times > 0) & (loss_times >= tcut) & (loss_times < trace_time)
    confined = (loss_times < 0) | (loss_times >= trace_time)

    # Background: regular orbits (green)
    ax.scatter(theta[regular & confined], pitch[regular & confined],
               c="tab:green", s=4, marker="s", alpha=0.6, edgecolors="none",
               label="Regular (confined)")

    # Regular but lost (false negatives -- should be rare)
    regular_lost = regular & ~confined
    if np.sum(regular_lost) > 0:
        ax.scatter(theta[regular_lost], pitch[regular_lost],
                   c="tab:green", s=4, marker="s", alpha=0.3, edgecolors="none")

    # Early losses (blue circles)
    early = prompt_loss | (chaotic & lost_early)
    ax.scatter(theta[early], pitch[early],
               c="tab:blue", s=8, marker="o", alpha=0.7, edgecolors="none",
               label="Early loss ($t < t_\\mathrm{cut}$)")

    # Chaotic late losses (red crosses)
    late = chaotic & lost_late
    ax.scatter(theta[late], pitch[late],
               c="tab:red", s=12, marker="x", alpha=0.7, linewidths=0.8,
               label="Chaotic (late loss)")

    # Chaotic but confined (yellow -- false positives)
    chaotic_conf = chaotic & confined
    if np.sum(chaotic_conf) > 0:
        ax.scatter(theta[chaotic_conf], pitch[chaotic_conf],
                   c="gold", s=8, marker="x", alpha=0.5, linewidths=0.8,
                   label="Chaotic (confined)")

    # Trapped-passing boundary: theta_trap = 0 means pitch^2 = 1 - B/B_max
    # Mark approximate boundary from trap_par sign change
    ax.axhline(y=0, color="white", linewidth=1.5, alpha=0.8)

    ax.set_xlabel(r"$\vartheta / \pi$")
    ax.set_ylabel(r"$v_\parallel / v$")
    ax.set_xlim([0, 2])
    ax.set_ylim([-1, 1])
    ax.legend(loc="upper right", fontsize=8, markerscale=1.5)
    ax.set_title("Orbit classification at $s = 0.6$")

    fig.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_volume_classification(data_dir, output_path):
    """
    Volume classification plot: (s, J_perp) colored by orbit type.

    Reuses the logic from plot_classification.py.
    """
    data_dir = Path(data_dir)

    class_data = np.loadtxt(data_dir / "class_parts.dat")
    s = class_data[:, 1]
    perp_inv = class_data[:, 2]
    icl_jpar = class_data[:, 3].astype(int)

    bminmax = np.loadtxt(data_dir / "bminmax.dat")

    ns, nperp = 50, 100
    hs = 1.0 / ns
    pmax = np.max(perp_inv)
    hp = pmax / nperp

    prompt = np.zeros((nperp, ns))
    regular = np.zeros((nperp, ns))
    stochastic = np.zeros((nperp, ns))

    for ipart in range(len(s)):
        i = min(ns, max(1, int(np.ceil(s[ipart] / hs)))) - 1
        k = min(nperp, max(1, int(np.ceil(perp_inv[ipart] / hp)))) - 1
        if icl_jpar[ipart] == 0:
            prompt[k, i] += 1.0
        elif icl_jpar[ipart] == 1:
            regular[k, i] += 1.0
        elif icl_jpar[ipart] == 2:
            stochastic[k, i] += 1.0

    tot = prompt + regular + stochastic
    cla = regular - prompt
    cla[tot == 0] = np.nan
    tot[tot == 0] = np.nan
    cla = cla / tot

    bmin_global = np.min(bminmax[:, 1])

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(cla, origin="lower", extent=[0, 1, 0, 1], aspect="auto",
                   vmin=-1.0, vmax=1.0)
    ax.plot(bminmax[:, 0], bmin_global / bminmax[:, 1], "k", lw=1.5)
    ax.plot(bminmax[:, 0], bmin_global / bminmax[:, 2], "k--", lw=1.0)

    jpmin = bmin_global / np.max(bminmax[:, 2])
    ax.set_ylim(jpmin, 1)
    ax.set_xlabel(r"Normalized toroidal flux $s$")
    ax.set_ylabel(r"Perpendicular invariant $J_\perp$")
    ax.set_title("Volume classification ($J_\\parallel$ classifier)")

    cb = plt.colorbar(im, ax=ax)
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels(["Prompt losses", "Chaotic", "Regular"])

    fig.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def print_table1(data_s03, data_s06, class_dir_s03, class_dir_s06):
    """
    Table 1 style: fractions of regular orbits in trapped and passing regions.
    """
    print("\n  Table 1: Fractions of regular orbits")
    print("  " + "-" * 55)
    print(f"  {'Surface':>10s}  {'Regular trapped':>16s}  {'Regular passing':>16s}")
    print("  " + "-" * 55)

    for label, data, class_dir in [("s = 0.3", data_s03, class_dir_s03),
                                    ("s = 0.6", data_s06, class_dir_s06)]:
        class_data = np.loadtxt(Path(class_dir) / "class_parts.dat")
        icl_jpar = class_data[:, 3].astype(int)

        trap_par = data.trap_parameter
        trapped = trap_par >= 0
        passing = trap_par < 0

        n_trapped = np.sum(trapped)
        n_passing = np.sum(passing)

        regular_trapped = np.sum((icl_jpar == 1) & trapped)
        regular_passing = np.sum((icl_jpar == 1) & passing)

        frac_trapped = regular_trapped / n_trapped if n_trapped > 0 else 0.0
        frac_passing = regular_passing / n_passing if n_passing > 0 else 0.0

        print(f"  {label:>10s}  {frac_trapped:>16.4f}  {frac_passing:>16.4f}")

    print("  " + "-" * 55)


def print_loss_summary(data, label):
    """Print loss summary for a single-surface run."""
    n = data.n_particles
    lost = data.lost_mask
    n_lost = np.sum(lost)
    fc = 1.0 - n_lost / n
    sigma = np.sqrt(fc * (1 - fc) / n)
    print(f"  {label}: {n} particles, {n_lost} lost, "
          f"f_c = {fc:.4f} +/- {1.96*sigma:.4f}")


def main():
    if len(sys.argv) > 1:
        base_dir = Path(sys.argv[1])
    else:
        base_dir = Path("/tmp/simple_classification")

    output_dir = base_dir / "plots"
    output_dir.mkdir(exist_ok=True)

    dir_s03 = base_dir / "s03"
    dir_s06 = base_dir / "s06"
    dir_vol = base_dir / "volume"

    # Check which runs are available
    have_s03 = (dir_s03 / "times_lost.dat").exists()
    have_s06 = (dir_s06 / "times_lost.dat").exists()
    have_vol = (dir_vol / "class_parts.dat").exists()

    if not (have_s03 or have_s06 or have_vol):
        print("No results found. Run run_all.sh first.")
        sys.exit(1)

    print("Loading data...")
    data_s03 = load_loss_data(dir_s03) if have_s03 else None
    data_s06 = load_loss_data(dir_s06) if have_s06 else None

    # Print loss summaries
    if data_s03 is not None:
        print_loss_summary(data_s03, "s=0.3")
    if data_s06 is not None:
        print_loss_summary(data_s06, "s=0.6")

    # Figure 6 style: loss time vs trapping parameter
    if data_s03 is not None and data_s06 is not None:
        print("\nGenerating Figure 6 style loss plots...")
        plot_paper_loss_figure(data_s03, data_s06,
                               output_dir / "fig6_losses_qh.pdf")

    # Individual KDE plots (reusing pysimple.plotting)
    if data_s03 is not None:
        print("Generating KDE loss density for s=0.3...")
        plot_kde_loss_density(data_s03,
                              output_path=output_dir / "kde_s03.pdf",
                              show=False)
        print(f"  Saved: {output_dir / 'kde_s03.pdf'}")

    if data_s06 is not None:
        print("Generating KDE loss density for s=0.6...")
        plot_kde_loss_density(data_s06,
                              output_path=output_dir / "kde_s06.pdf",
                              show=False)
        print(f"  Saved: {output_dir / 'kde_s06.pdf'}")

    # Figure 8 style: classification in (theta, pitch) space
    if have_s06:
        print("\nGenerating Figure 8 style classification plot...")
        plot_classification_theta_pitch(
            dir_s06, output_dir / "fig8_classification_qh.pdf")

    # Volume classification plot
    if have_vol:
        print("\nGenerating volume classification plot...")
        plot_volume_classification(dir_vol, output_dir / "volume_classification.pdf")

    # Table 1 style output
    if data_s03 is not None and data_s06 is not None:
        print_table1(data_s03, data_s06, dir_s03, dir_s06)

    print(f"\nAll plots saved to: {output_dir}")


if __name__ == "__main__":
    main()
