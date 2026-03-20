#!/usr/bin/env python3
"""
Reproduce all classification figures from the JPP 2020 paper:
  Albert, Kasilov, Kernbichler, "Accelerated methods for direct computation
  of fusion alpha particle losses within stellarator optimization",
  J. Plasma Phys. 86, 815860201 (2020).

Uses the actual paper equilibria:
  - QI: Drevlak et al. 2014 (wout_23_1900_fix_bdry.nc, nfp=5)
  - QH: Drevlak et al. 2018 (wout_qh_8_7.nc, nfp=5)
  - QA: Henneberg et al. 2019 (wout_henneberg_qa.nc, nfp=2)

Generates:
  - Figure 5: QI loss plots at s=0.3 and s=0.6
  - Figure 6: QH loss plots at s=0.3 and s=0.6
  - Figure 7: QA loss plots at s=0.3 and s=0.6
  - Figure 8: Orbit classification in (theta, v_par/v) for all 3 configs
  - Table 1: Regular orbit fractions for trapped/passing regions
  - Volume classification plots for QH, QI, QA

Usage:
    python plot_paper_results.py [base_dir]

    base_dir defaults to /tmp/simple_classification
"""

import sys
from pathlib import Path

import numpy as np

# Import pysimple.plotting without requiring compiled backend
_plotting_path = (
    Path(__file__).resolve().parent.parent.parent / "python" / "pysimple" / "plotting.py"
)
import importlib.util

_spec = importlib.util.spec_from_file_location("pysimple.plotting", _plotting_path)
_plotting = importlib.util.module_from_spec(_spec)
sys.modules["pysimple.plotting"] = _plotting
_spec.loader.exec_module(_plotting)

load_loss_data = _plotting.load_loss_data

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Configuration labels and directory name mappings
CONFIGS = {
    "QI": {"label": "QI (Drevlak 2014)", "s03": "qi_s03", "s06": "qi_s06", "fig8": "qi_fig8"},
    "QH": {"label": "QH (Drevlak 2018)", "s03": "qh_s03", "s06": "qh_s06", "fig8": "qh_fig8"},
    "QA": {"label": "QA (Henneberg 2019)", "s03": "qa_s03", "s06": "qa_s06", "fig8": "qa_fig8"},
}
FIGURE_NUMBERS = {"QI": 5, "QH": 6, "QA": 7}


def plot_loss_panel(ax, data, title):
    """Plot a single loss panel: KDE density + dots + confined fraction.

    Matches the style of Figures 5-7 in the paper: black dots for lost
    particles in (log10(t), trapping parameter) space with blue KDE shading,
    a red confined fraction curve on the right axis, and error bands.
    """
    lost = data.lost_mask
    tlost = np.abs(data.loss_times[lost])
    trap_par = data.trap_parameter[lost]
    t_threshold = 0.99 * data.trace_time
    kde_mask = tlost < t_threshold

    # KDE density shading
    if np.sum(kde_mask) > 10:
        log_tlost = np.log10(tlost[kde_mask])
        trap_kde = trap_par[kde_mask]
        try:
            kde = gaussian_kde([log_tlost, trap_kde])
            xg, yg = np.mgrid[-5:0:100j, -2:1:150j]
            positions = np.vstack([xg.ravel(), yg.ravel()])
            density = np.reshape(kde(positions).T, xg.shape)
            ax.imshow(
                np.rot90(density), cmap=plt.cm.Blues,
                extent=[-5, 0, -2, 1], aspect="auto",
            )
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
         r"$10^{-2}$", r"$10^{-1}$", r"$1$"]
    )
    ax.set_title(title)

    # Confined fraction on twin axis
    ax2 = ax.twinx()
    n = data.n_particles
    t_grid = np.logspace(-5, 0, 200)
    fc = np.array([np.sum(np.abs(data.loss_times) >= t) / n for t in t_grid])
    sigma = np.sqrt(fc * (1 - fc) / n)

    ax2.plot(np.log10(t_grid), fc, color="tab:red", lw=1.5)
    ax2.fill_between(
        np.log10(t_grid), fc - 1.96 * sigma, fc + 1.96 * sigma,
        color="tab:red", alpha=0.15,
    )
    ax2.set_ylim([0.6, 1.05])
    ax2.set_ylabel(r"confined fraction $f_c$", color="tab:red")
    ax2.tick_params("y", colors="tab:red")


def plot_loss_figure(data_s03, data_s06, config_key, output_path):
    """Generate a two-panel loss figure (Figures 5-7 style)."""
    fig_num = FIGURE_NUMBERS[config_key]
    label = CONFIGS[config_key]["label"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    plot_loss_panel(axes[0], data_s03, f"({chr(97)}) $s = 0.3$")
    plot_loss_panel(axes[1], data_s06, f"({chr(98)}) $s = 0.6$")

    fig.suptitle(f"Figure {fig_num}: {label}", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def _read_namelist_param(simple_in, param_name, default):
    """Read a single parameter value from a simple.in namelist file."""
    if not simple_in.exists():
        return default
    for line in simple_in.read_text().splitlines():
        if param_name in line and "=" in line:
            if param_name == "tcut" and "ntcut" in line:
                continue
            val = line.split("=")[1].split("!")[0].strip().rstrip(",")
            return float(val.replace("d", "e"))
    return default


def _classify_particles(fig8_dir):
    """Load and classify particles from a Fig 8 run directory.

    Returns (theta, pitch, category) where category is:
      0 = regular (confined), 1 = early loss, 2 = chaotic/late loss
    Returns None if data files are missing.
    """
    for required in ("class_parts.dat", "times_lost.dat", "start.dat"):
        if not (fig8_dir / required).exists():
            return None

    class_data = np.loadtxt(fig8_dir / "class_parts.dat")
    iclass_jpar = class_data[:, 3].astype(int)

    tl_data = np.loadtxt(fig8_dir / "times_lost.dat")
    loss_times = tl_data[:, 1]

    start_data = np.loadtxt(fig8_dir / "start.dat")
    theta = np.mod(start_data[:, 1], 2 * np.pi) / np.pi
    pitch = start_data[:, 4]

    simple_in = fig8_dir / "simple.in"
    trace_time = _read_namelist_param(simple_in, "trace_time", 1.0)
    tcut = _read_namelist_param(simple_in, "tcut", 0.1)

    regular = iclass_jpar == 1
    prompt_loss = iclass_jpar == 0
    chaotic = iclass_jpar == 2

    lost_early = (loss_times > 0) & (loss_times < tcut)
    lost_late = (loss_times > 0) & (loss_times >= tcut) & (loss_times < trace_time)
    confined = (loss_times < 0) | (loss_times >= trace_time)

    # Category: 0=regular, 1=early loss, 2=chaotic
    category = np.full(len(theta), -1, dtype=int)
    category[regular] = 0
    category[prompt_loss | (chaotic & lost_early)] = 1
    category[chaotic & (lost_late | confined)] = 2
    # Regular but lost -- still mark as regular for background
    category[regular & ~confined] = 0

    return theta, pitch, category


def _plot_fig8_panel(ax, fig8_dir, panel_label):
    """Plot a single Fig 8 panel with dense imshow background and overlay markers.

    Uses a coarse grid with Gaussian smoothing so the background fills the
    entire (theta/pi, v_par/v) plane -- the "brazilian flag" look from the
    paper.  Colors: green=regular, blue=early loss, yellow=chaotic.
    Overlay markers: circles for early losses, crosses for chaotic.
    """
    from scipy.ndimage import gaussian_filter

    result = _classify_particles(fig8_dir)
    if result is None:
        ax.set_title(f"{panel_label} -- no data")
        return

    theta, pitch, category = result

    n_theta, n_pitch = 50, 50
    theta_edges = np.linspace(0, 2, n_theta + 1)
    pitch_edges = np.linspace(-1, 1, n_pitch + 1)

    # Count particles of each type in each grid cell
    counts = np.zeros((3, n_pitch, n_theta))
    for cat_id in range(3):
        mask = category == cat_id
        if np.any(mask):
            h, _, _ = np.histogram2d(
                pitch[mask], theta[mask],
                bins=[pitch_edges, theta_edges],
            )
            counts[cat_id] = h

    # Smooth counts to fill gaps and produce continuous colored regions
    sigma = 1.5
    smoothed = np.array([gaussian_filter(counts[i], sigma) for i in range(3)])
    total_smooth = smoothed.sum(axis=0)

    # Build RGB image from smoothed fractional composition
    color_map = np.array([
        [0.2, 0.7, 0.2],   # green (regular)
        [0.2, 0.4, 0.9],   # blue (early loss)
        [0.9, 0.8, 0.1],   # yellow (chaotic)
    ])
    rgb = np.ones((n_pitch, n_theta, 3))
    filled = total_smooth > 1e-10
    for cat_id in range(3):
        frac = np.where(filled, smoothed[cat_id] / np.where(filled, total_smooth, 1.0), 0.0)
        for ch in range(3):
            rgb[:, :, ch] = np.where(
                filled,
                rgb[:, :, ch] - frac * (1.0 - color_map[cat_id, ch]),
                rgb[:, :, ch],
            )
    rgb = np.clip(rgb, 0, 1)

    ax.imshow(
        rgb, origin="lower", extent=[0, 2, -1, 1], aspect="auto",
        interpolation="bilinear",
    )

    # Overlay markers for early losses (circles) and chaotic (crosses)
    early = category == 1
    if np.any(early):
        ax.scatter(
            theta[early], pitch[early],
            c="blue", s=2, marker="o", alpha=0.3, edgecolors="none",
            label=r"Early loss ($t < t_\mathrm{cut}$)",
        )

    chaotic_mask = category == 2
    if np.any(chaotic_mask):
        ax.scatter(
            theta[chaotic_mask], pitch[chaotic_mask],
            c="red", s=4, marker="x", alpha=0.4, linewidths=0.5,
            label="Chaotic",
        )

    # Trapped-passing boundary
    ax.axhline(y=0, color="white", linewidth=1.5, alpha=0.8)

    ax.set_xlabel(r"$\vartheta / \pi$")
    ax.set_xlim([0, 2])
    ax.set_ylim([-1, 1])
    ax.set_title(panel_label)


def plot_fig8(base_dir, output_path):
    """Generate Figure 8: orbit classification in (theta, v_par/v) space.

    Three panels (a, b, c) for QI, QH, QA at s=0.6, phi = phi_n/2.
    Uses a dense colored grid (imshow) for the background with overlay markers.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    panel_labels = ["(a) QI", "(b) QH", "(c) QA"]

    for ax, config_key, panel_label in zip(axes, ["QI", "QH", "QA"], panel_labels):
        fig8_dir = base_dir / CONFIGS[config_key]["fig8"]
        _plot_fig8_panel(ax, fig8_dir, panel_label)

    # Only the leftmost panel gets a y-label
    axes[0].set_ylabel(r"$v_\parallel / v$")

    # Shared legend from whichever panel has handles
    handles, labels = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        for hi, li in zip(h, l):
            if li not in labels:
                handles.append(hi)
                labels.append(li)
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8,
                   bbox_to_anchor=(0.5, -0.05))

    fig.suptitle("Figure 8: Orbit classification at $s = 0.6$", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_volume_classification(data_dir, output_path, config_label=""):
    """Volume classification plot: (s, J_perp) colored by orbit type.

    Bins trapped particles only (iclass != 1 or trap_par >= 0) into a
    (s, J_perp) grid. Passing particles are shown as white background.
    The trapped-passing and deeply-trapped boundaries from bminmax.dat
    are overlaid.
    """
    data_dir = Path(data_dir)

    class_data = np.loadtxt(data_dir / "class_parts.dat")
    s = class_data[:, 1]
    perp_inv = class_data[:, 2]
    icl_jpar = class_data[:, 3].astype(int)

    tl_data = np.loadtxt(data_dir / "times_lost.dat")
    trap_par = tl_data[:, 2]

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
        is_trapped = trap_par[ipart] >= 0
        if not is_trapped:
            continue
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
    im = ax.imshow(
        cla, origin="lower", extent=[0, 1, 0, 1], aspect="auto",
        vmin=-1.0, vmax=1.0,
    )
    ax.plot(bminmax[:, 0], bmin_global / bminmax[:, 1], "k", lw=1.5)
    ax.plot(bminmax[:, 0], bmin_global / bminmax[:, 2], "k--", lw=1.0)

    jpmin = bmin_global / np.max(bminmax[:, 2])
    ax.set_ylim(jpmin, 1)
    ax.set_xlabel(r"Normalized toroidal flux $s$")
    ax.set_ylabel(r"Perpendicular invariant $J_\perp$")
    title_suffix = f" -- {config_label}" if config_label else ""
    ax.set_title(f"Volume classification{title_suffix}")

    cb = plt.colorbar(im, ax=ax)
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels(["Prompt losses", "Chaotic", "Regular"])

    fig.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def print_table1(base_dir):
    """Print Table 1: fractions of regular orbits in trapped/passing regions."""
    print("\n  Table 1: Fractions of regular orbits")
    print("  " + "-" * 60)
    print(f"  {'Type':>6s}  {'Surface s':>10s}  {'Regular trapped':>16s}  {'Regular passing':>16s}")
    print("  " + "-" * 60)

    for config_key in ["QI", "QH", "QA"]:
        for surface, dir_name in [("0.3", CONFIGS[config_key]["s03"]),
                                   ("0.6", CONFIGS[config_key]["s06"])]:
            run_dir = base_dir / dir_name
            class_file = run_dir / "class_parts.dat"
            tl_file = run_dir / "times_lost.dat"

            if not class_file.exists() or not tl_file.exists():
                print(f"  {config_key:>6s}  {surface:>10s}  {'N/A':>16s}  {'N/A':>16s}")
                continue

            class_data = np.loadtxt(class_file)
            icl_jpar = class_data[:, 3].astype(int)

            tl_data = np.loadtxt(tl_file)
            trap_par = tl_data[:, 2]

            trapped = trap_par >= 0
            passing = trap_par < 0

            n_trapped = np.sum(trapped)
            n_passing = np.sum(passing)

            regular_trapped = np.sum((icl_jpar == 1) & trapped)
            regular_passing = np.sum((icl_jpar == 1) & passing)

            frac_trapped = regular_trapped / n_trapped if n_trapped > 0 else 0.0
            frac_passing = regular_passing / n_passing if n_passing > 0 else 0.0

            print(f"  {config_key:>6s}  {surface:>10s}  {frac_trapped:>16.4f}  {frac_passing:>16.4f}")

    print("  " + "-" * 60)


def print_loss_summary(data, label):
    """Print loss summary for a single-surface run."""
    n = data.n_particles
    lost = data.lost_mask
    n_lost = np.sum(lost)
    fc = 1.0 - n_lost / n
    sigma = np.sqrt(fc * (1 - fc) / n)
    print(f"  {label}: {n} particles, {n_lost} lost, "
          f"f_c = {fc:.4f} +/- {1.96 * sigma:.4f}")


def main():
    if len(sys.argv) > 1:
        base_dir = Path(sys.argv[1])
    else:
        base_dir = Path("/tmp/simple_classification")

    output_dir = base_dir / "plots"
    output_dir.mkdir(exist_ok=True)

    # Load data for all configurations
    datasets = {}
    for config_key in ["QI", "QH", "QA"]:
        for surface in ["s03", "s06"]:
            dir_name = CONFIGS[config_key][surface]
            run_dir = base_dir / dir_name
            tl_file = run_dir / "times_lost.dat"
            key = f"{config_key}_{surface}"
            if tl_file.exists():
                datasets[key] = load_loss_data(run_dir)
                s_label = "0.3" if surface == "s03" else "0.6"
                print_loss_summary(datasets[key], f"{config_key} s={s_label}")

    # Figures 5-7: loss plots for each configuration
    for config_key, fig_num in FIGURE_NUMBERS.items():
        key_s03 = f"{config_key}_s03"
        key_s06 = f"{config_key}_s06"
        if key_s03 in datasets and key_s06 in datasets:
            print(f"\nGenerating Figure {fig_num} ({config_key} loss plots)...")
            plot_loss_figure(
                datasets[key_s03], datasets[key_s06], config_key,
                output_dir / f"fig{fig_num}_losses_{config_key.lower()}.pdf",
            )

    # Figure 8: orbit classification
    have_any_fig8 = any(
        (base_dir / CONFIGS[ck]["fig8"] / "class_parts.dat").exists()
        for ck in ["QI", "QH", "QA"]
    )
    if have_any_fig8:
        print("\nGenerating Figure 8 (orbit classification)...")
        plot_fig8(base_dir, output_dir / "fig8_classification.pdf")

    # Volume classification for all 3 configs
    volume_configs = [
        ("QH", "qh_volume"),
        ("QI", "qi_volume"),
        ("QA", "qa_volume"),
    ]
    for vol_label, vol_name in volume_configs:
        vol_dir = base_dir / vol_name
        if (vol_dir / "class_parts.dat").exists():
            out_name = f"volume_classification_{vol_label.lower()}.pdf"
            print(f"\nGenerating volume classification plot ({vol_label})...")
            plot_volume_classification(vol_dir, output_dir / out_name, vol_label)

    # Table 1
    print_table1(base_dir)

    print(f"\nAll plots saved to: {output_dir}")


if __name__ == "__main__":
    main()
