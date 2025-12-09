"""
Loss statistics visualization for SIMPLE particle tracing.

This module provides publication-quality plots for analyzing particle losses
in stellarator simulations, including:
- Confined fraction over time (particles and energy)
- KDE density plots of loss time vs trapping parameter
- Energy loss distribution vs perpendicular invariant J_perp
- Starting position visualization

Based on visualize3-MK.py by Christopher Albert, Majid Khan, and Daniele Corrias.

Data File Formats
-----------------
The plotting functions read standard SIMPLE output files:

**start.dat** - Initial particle conditions
    Column 1: r (normalized toroidal flux s)
    Column 2: theta_vmec (poloidal angle)
    Column 3: varphi_vmec (toroidal angle)
    Column 4: normalized velocity module v/v0 (=1 for monoenergetic)
    Column 5: pitch angle cosine v_parallel/v (lambda)

**times_lost.dat** - Loss times and final states
    Column 1: Particle index (corresponds to line in start.dat)
    Column 2: Loss time t_loss [s]
              - trace_time if confined/regular
              - -1 if skipped due to contr_pp (deep passing)
    Column 3: Trapping parameter eta
              - eta=1: deeply trapped
              - eta=0: trapped-passing boundary
              - eta<0: passing particles
              See Eq. (3.1) in Accelerated Methods paper.
    Column 4: Starting s (normalized toroidal flux)
    Column 5: J_perp (perpendicular adiabatic invariant)
    Columns 6-10: Final state zend(1:5)
              - zend(1): final s
              - zend(2): final theta
              - zend(3): final phi
              - zend(4): final p = v/v0 (normalized velocity)
              - zend(5): final v_parallel/v

**confined_fraction.dat** - Time evolution of confinement
    Column 1: Time [s]
    Column 2: confpart_pass (confined passing fraction)
    Column 3: confpart_trap (confined trapped fraction)
    Column 4: Total number of particles
    Note: Total confined fraction = Column 2 + Column 3

Example
-------
>>> from pysimple.plotting import load_loss_data, plot_loss_statistics
>>> data = load_loss_data("./run_output/")
>>> fig = plot_loss_statistics(data, output_path="losses.png")
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from scipy.stats import gaussian_kde

    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False
    Figure = None


@dataclass
class LossData:
    """
    Container for SIMPLE loss statistics data.

    Attributes
    ----------
    n_particles : int
        Total number of particles
    loss_times : np.ndarray
        Loss time for each particle [s]. Equal to trace_time if confined.
    trap_parameter : np.ndarray
        Trapping parameter eta. Positive=trapped, negative=passing.
    perp_invariant : np.ndarray
        Perpendicular adiabatic invariant J_perp.
    start_s : np.ndarray
        Starting normalized toroidal flux.
    start_theta : np.ndarray
        Starting poloidal angle [rad].
    start_phi : np.ndarray
        Starting toroidal angle [rad].
    start_pitch : np.ndarray
        Starting pitch angle cosine v_parallel/v.
    final_p : np.ndarray
        Final normalized velocity p = v/v0 at loss time.
        Energy at loss = final_p**2 (normalized to birth energy).
    time_grid : np.ndarray
        Time grid from confined_fraction.dat [s].
    confined_pass : np.ndarray
        Confined passing fraction at each time.
    confined_trap : np.ndarray
        Confined trapped fraction at each time.
    trace_time : float
        Maximum tracing time [s].
    """

    n_particles: int
    loss_times: np.ndarray
    trap_parameter: np.ndarray
    perp_invariant: np.ndarray
    start_s: np.ndarray
    start_theta: np.ndarray
    start_phi: np.ndarray
    start_pitch: np.ndarray
    final_p: np.ndarray
    time_grid: np.ndarray
    confined_pass: np.ndarray
    confined_trap: np.ndarray
    trace_time: float

    @property
    def confined_fraction(self) -> np.ndarray:
        """Total confined fraction (passing + trapped) at each time."""
        return self.confined_pass + self.confined_trap

    @property
    def lost_mask(self) -> np.ndarray:
        """Boolean mask for particles that were lost (not confined)."""
        return (self.loss_times > 0) & (self.loss_times < self.trace_time)

    @property
    def confined_mask(self) -> np.ndarray:
        """Boolean mask for particles that remained confined."""
        return self.loss_times >= self.trace_time

    @property
    def skipped_mask(self) -> np.ndarray:
        """Boolean mask for particles skipped (deep passing, contr_pp)."""
        return self.loss_times < 0


def load_loss_data(directory: str | Path) -> LossData:
    """
    Load SIMPLE loss statistics from output directory.

    Parameters
    ----------
    directory : str or Path
        Directory containing times_lost.dat, confined_fraction.dat, and start.dat

    Returns
    -------
    LossData
        Container with all loss statistics data.

    Raises
    ------
    FileNotFoundError
        If required files are not found in directory.
    """
    directory = Path(directory)

    # Load times_lost.dat
    # Columns: index, t_lost, trap_par, s_start, J_perp, zend(1:5)
    times_lost_file = directory / "times_lost.dat"
    if not times_lost_file.exists():
        raise FileNotFoundError(f"times_lost.dat not found in {directory}")

    tl_data = np.loadtxt(times_lost_file)
    n_particles = tl_data.shape[0]

    loss_times = tl_data[:, 1]
    trap_parameter = tl_data[:, 2]
    start_s_tl = tl_data[:, 3]
    perp_invariant = tl_data[:, 4]

    # Final state: zend columns 5-9 (0-indexed: columns 5,6,7,8,9)
    # zend(4) = final normalized velocity p = v/v0
    if tl_data.shape[1] >= 10:
        final_p = tl_data[:, 8]  # zend(4) is column index 8
    else:
        # Old format without zend - assume p=1 (no collisions)
        final_p = np.ones(n_particles)

    # Load start.dat for initial conditions
    start_file = directory / "start.dat"
    if start_file.exists():
        start_data = np.loadtxt(start_file)
        # Handle case where start.dat may have more particles than times_lost.dat
        start_data = start_data[:n_particles]
        start_s = start_data[:, 0]
        start_theta = np.mod(start_data[:, 1], 2 * np.pi)
        start_phi = np.mod(start_data[:, 2], 2 * np.pi)
        start_pitch = start_data[:, 4]
    else:
        # Use s from times_lost.dat, set others to zero
        start_s = start_s_tl
        start_theta = np.zeros(n_particles)
        start_phi = np.zeros(n_particles)
        start_pitch = np.zeros(n_particles)

    # Load confined_fraction.dat
    # Columns: time, confpart_pass, confpart_trap, ntestpart
    conf_file = directory / "confined_fraction.dat"
    if conf_file.exists():
        conf_data = np.loadtxt(conf_file)
        # Skip first row (t=0, all confined)
        time_grid = np.abs(conf_data[1:, 0])
        confined_pass = conf_data[1:, 1]
        confined_trap = conf_data[1:, 2]
    else:
        # Generate from loss times if file not available
        time_grid = np.logspace(-5, 0, 100)
        confined_pass = np.zeros_like(time_grid)
        confined_trap = np.zeros_like(time_grid)
        for i, t in enumerate(time_grid):
            confined = np.abs(loss_times) >= t
            passing = trap_parameter < 0
            confined_pass[i] = np.sum(confined & passing) / n_particles
            confined_trap[i] = np.sum(confined & ~passing) / n_particles

    # Determine trace_time from maximum loss time or last time in grid
    trace_time = max(np.max(np.abs(loss_times)), time_grid[-1] if len(time_grid) > 0 else 1.0)

    return LossData(
        n_particles=n_particles,
        loss_times=loss_times,
        trap_parameter=trap_parameter,
        perp_invariant=perp_invariant,
        start_s=start_s,
        start_theta=start_theta,
        start_phi=start_phi,
        start_pitch=start_pitch,
        final_p=final_p,
        time_grid=time_grid,
        confined_pass=confined_pass,
        confined_trap=confined_trap,
        trace_time=trace_time,
    )


def compute_energy_confined_fraction(
    data: LossData,
    time_grid: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute energy confined fraction over time.

    The energy confined fraction accounts for the fact that particles
    lost early carry away more energy than those lost late (after
    slowing down via collisions).

    Energy at loss time = final_p**2 where final_p = v/v0 is the
    normalized velocity from times_lost.dat column 9.

    Parameters
    ----------
    data : LossData
        Loss statistics data from load_loss_data().
    time_grid : np.ndarray, optional
        Time points to evaluate. Defaults to data.time_grid.

    Returns
    -------
    times : np.ndarray
        Time grid [s].
    energy_fraction : np.ndarray
        Energy confined fraction at each time (0 to 1).
    """
    if time_grid is None:
        time_grid = data.time_grid

    # Total initial energy (normalized: each particle has E=1 at birth)
    e_total = float(data.n_particles)

    # Energy carried away by each lost particle = final_p**2
    # For confined particles or skipped particles, they don't contribute to loss
    energy_at_loss = data.final_p**2

    energy_fraction = np.zeros_like(time_grid)

    for i, t in enumerate(time_grid):
        # Particles lost before time t
        lost_before_t = (data.loss_times > 0) & (data.loss_times < t)

        # Energy lost = sum of energy carried away by lost particles
        e_lost = np.sum(energy_at_loss[lost_before_t])

        # Energy confined = total - lost
        energy_fraction[i] = (e_total - e_lost) / e_total

    return time_grid, energy_fraction


def compute_energy_loss_distribution(
    data: LossData,
    n_bins: int = 50,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute energy loss distribution binned by perpendicular invariant J_perp.

    Parameters
    ----------
    data : LossData
        Loss statistics data.
    n_bins : int
        Number of J_perp bins.

    Returns
    -------
    jperp_bins : np.ndarray
        J_perp bin centers (normalized to max).
    particle_count : np.ndarray
        Number of particles in each bin.
    energy_lost : np.ndarray
        Total energy lost in each bin (sum of final_p**2 for lost particles).
    """
    lost = data.lost_mask
    jperp_lost = data.perp_invariant[lost]
    energy_lost_per_particle = data.final_p[lost] ** 2

    if len(jperp_lost) == 0:
        return np.zeros(n_bins), np.zeros(n_bins), np.zeros(n_bins)

    # Normalize J_perp to [0, 1]
    jperp_max = np.max(data.perp_invariant)
    if jperp_max > 0:
        jperp_norm = jperp_lost / jperp_max
    else:
        jperp_norm = jperp_lost

    # Bin edges
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Histogram
    particle_count, _ = np.histogram(jperp_norm, bins=bin_edges)
    energy_lost, _ = np.histogram(jperp_norm, bins=bin_edges, weights=energy_lost_per_particle)

    return bin_centers, particle_count.astype(float), energy_lost


def plot_loss_statistics(
    data: LossData,
    output_path: Optional[str | Path] = None,
    show: bool = True,
    figsize: Tuple[float, float] = (10, 8),
    dpi: int = 150,
) -> Figure:
    """
    Create comprehensive loss statistics figure with 4 panels.

    Panel layout:
    +----------------------------------+----------------------------------+
    | (a) Confined Fraction vs Time    | (b) KDE: Loss Time vs Trap Param |
    |     - Red: particle fraction     |     - Blue KDE density           |
    |     - Blue: energy fraction      |     - Red: confined fraction     |
    +----------------------------------+----------------------------------+
    | (c) Energy Loss vs J_perp        | (d) Starting Positions           |
    |                                  |     θ vs λ colored by loss time  |
    +----------------------------------+----------------------------------+

    Parameters
    ----------
    data : LossData
        Loss statistics data from load_loss_data().
    output_path : str or Path, optional
        If provided, save figure to this path.
    show : bool
        If True, display the figure.
    figsize : tuple
        Figure size in inches (width, height).
    dpi : int
        Resolution for saved figure.

    Returns
    -------
    Figure
        Matplotlib figure object.

    Raises
    ------
    ImportError
        If matplotlib or scipy are not available.
    """
    if not HAS_PLOTTING:
        raise ImportError(
            "Plotting requires matplotlib and scipy. "
            "Install with: pip install matplotlib scipy"
        )

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # =========================================================================
    # Panel (a): Confined fraction vs time - particles AND energy
    # =========================================================================
    ax_conf = axes[0, 0]

    # Particle confined fraction (from confined_fraction.dat)
    t_conf = data.time_grid
    conf_frac = data.confined_fraction

    # Energy confined fraction (computed from final velocities)
    t_energy, energy_frac = compute_energy_confined_fraction(data)

    # Plot with log scale on x-axis
    valid_t = t_conf > 0
    ax_conf.semilogx(t_conf[valid_t], conf_frac[valid_t], "r-", linewidth=1.5, label="Particles")
    ax_conf.semilogx(t_energy[valid_t], energy_frac[valid_t], "b-", linewidth=1.5, label="Energy")

    ax_conf.set_xlabel(r"Time $t$ / s")
    ax_conf.set_ylabel("Confined fraction")
    ax_conf.set_ylim([0, 1.05])
    ax_conf.legend(loc="lower left")
    ax_conf.grid(True, alpha=0.3)
    ax_conf.set_title("(a) Confined fraction vs time")

    # =========================================================================
    # Panel (b): KDE density plot - loss time vs trapping parameter
    # This is the "Brazilian flag" style plot from ISHW papers
    # =========================================================================
    ax_kde = axes[0, 1]

    # Filter to lost particles only (exclude confined and skipped)
    lost = data.lost_mask
    tlost = np.abs(data.loss_times[lost])
    trap_par = data.trap_parameter[lost]

    # Further filter to reasonable loss times for KDE
    t_threshold = 0.99 * data.trace_time
    kde_mask = tlost < t_threshold

    if np.sum(kde_mask) > 10:
        # Compute 2D KDE (Kernel Density Estimate)
        # This shows the density of particles in (log10(t_lost), trap_par) space
        log_tlost = np.log10(tlost[kde_mask])
        trap_kde = trap_par[kde_mask]

        try:
            kde = gaussian_kde([log_tlost, trap_kde])

            # Create evaluation grid
            # X: log10(time) from 10^-5 to 10^0 (1 second)
            # Y: trapping parameter from -2 to 1
            X, Y = np.mgrid[-5:0:100j, -2:1:150j]
            positions = np.vstack([X.ravel(), Y.ravel()])
            Z = np.reshape(kde(positions).T, X.shape)

            # Plot KDE as image
            ax_kde.imshow(
                np.rot90(Z),
                cmap=plt.cm.Blues,
                extent=[-5, 0, -2, 1],
                aspect="auto",
            )
        except np.linalg.LinAlgError:
            # KDE can fail if data is degenerate
            pass

    # Overlay scatter of individual particles
    if np.sum(lost) > 0:
        ax_kde.plot(
            np.log10(tlost),
            trap_par,
            "k,",
            markersize=0.5,
            alpha=0.5,
        )

    # Dashed line at trap_par = 0 (trapped-passing boundary)
    ax_kde.axhline(y=0, color="k", linestyle="--", linewidth=0.8)

    ax_kde.set_ylim([-2.1, 1.0])
    ax_kde.set_xlim([-5.0, 0.1])
    ax_kde.set_ylabel("Trapping parameter $\\eta$")
    ax_kde.set_xlabel(r"Loss time $t$ / s")

    # Custom x-axis labels showing powers of 10
    ax_kde.set_xticks([-5, -4, -3, -2, -1, 0])
    ax_kde.set_xticklabels(
        [r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$1$"]
    )

    # Add confined fraction on twin y-axis (red curve)
    ax_conf2 = ax_kde.twinx()
    valid_t = t_conf > 0
    ax_conf2.plot(np.log10(t_conf[valid_t]), conf_frac[valid_t], color="tab:red", linewidth=1.5)
    ax_conf2.set_ylim([0.5, 1.25])
    ax_conf2.set_ylabel("Confined fraction", color="tab:red")
    ax_conf2.tick_params("y", colors="tab:red")

    ax_kde.set_title("(b) Loss time vs trapping parameter (KDE)")

    # =========================================================================
    # Panel (c): Energy loss distribution vs J_perp
    # =========================================================================
    ax_energy = axes[1, 0]

    jperp_bins, particle_count, energy_lost = compute_energy_loss_distribution(data)

    # Normalize energy lost to get fraction
    total_particles = data.n_particles
    if total_particles > 0:
        energy_frac_per_bin = energy_lost / total_particles
        particle_frac_per_bin = particle_count / total_particles
    else:
        energy_frac_per_bin = energy_lost
        particle_frac_per_bin = particle_count

    ax_energy.bar(
        jperp_bins,
        energy_frac_per_bin,
        width=1.0 / len(jperp_bins),
        alpha=0.7,
        label="Energy lost",
        color="tab:blue",
    )
    ax_energy.bar(
        jperp_bins,
        particle_frac_per_bin,
        width=1.0 / len(jperp_bins),
        alpha=0.4,
        label="Particles lost",
        color="tab:red",
    )

    ax_energy.set_xlabel(r"$J_\perp$ (normalized)")
    ax_energy.set_ylabel("Lost fraction per bin")
    ax_energy.legend()
    ax_energy.set_title("(c) Loss distribution vs $J_\\perp$")

    # =========================================================================
    # Panel (d): Starting positions - theta vs pitch colored by loss time
    # =========================================================================
    ax_start = axes[1, 1]

    # Color by log10 of loss time
    # Confined particles (tlost = trace_time) will be yellow
    # Early losses will be dark (purple/blue)
    tlost_all = np.abs(data.loss_times)
    tlost_all = np.maximum(tlost_all, 1e-10)  # Avoid log(0)

    scatter = ax_start.scatter(
        data.start_theta,
        data.start_pitch,
        c=np.log10(tlost_all),
        s=1,
        cmap="viridis",
        vmin=-5,
        vmax=0,
    )

    cbar = plt.colorbar(scatter, ax=ax_start)
    cbar.set_label(r"$\log_{10}(t_{\mathrm{loss}})$")

    ax_start.set_xlabel(r"$\theta$ (poloidal angle)")
    ax_start.set_ylabel(r"$\lambda = v_\parallel / v$ (pitch)")
    ax_start.set_xlim([0, 2 * np.pi])
    ax_start.set_ylim([-1, 1])
    ax_start.set_title("(d) Starting positions colored by loss time")

    # =========================================================================
    # Finalize
    # =========================================================================
    plt.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_confined_fraction(
    data: LossData,
    output_path: Optional[str | Path] = None,
    show: bool = True,
    include_energy: bool = True,
) -> Figure:
    """
    Simple plot of confined fraction vs time.

    Parameters
    ----------
    data : LossData
        Loss statistics data.
    output_path : str or Path, optional
        If provided, save figure to this path.
    show : bool
        If True, display the figure.
    include_energy : bool
        If True, also plot energy confined fraction.

    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    if not HAS_PLOTTING:
        raise ImportError("Plotting requires matplotlib")

    fig, ax = plt.subplots(figsize=(6, 4))

    t_conf = data.time_grid
    valid_t = t_conf > 0

    ax.semilogx(
        t_conf[valid_t],
        data.confined_fraction[valid_t],
        "r-",
        linewidth=1.5,
        label="Particles",
    )

    if include_energy:
        t_energy, energy_frac = compute_energy_confined_fraction(data)
        ax.semilogx(
            t_energy[valid_t],
            energy_frac[valid_t],
            "b-",
            linewidth=1.5,
            label="Energy",
        )

    ax.set_xlabel(r"Time $t$ / s")
    ax.set_ylabel("Confined fraction")
    ax.set_ylim([0, 1.05])
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_energy_loss_vs_jperp(
    data_coll: LossData,
    data_nocoll: Optional[LossData] = None,
    output_path: Optional[str | Path] = None,
    show: bool = True,
    figsize: Tuple[float, float] = (7, 6),
    title: Optional[str] = None,
    xlim: Optional[float] = None,
    nperp: int = 100,
    smooth_sigma: float = 2.0,
) -> Figure:
    """
    Plot lost energy fraction vs J_perp in ISHW 2022 style.

    Uses the NEW approach: energy carried away by lost particles = final_p^2
    (where final_p = v/v0 at loss time from times_lost.dat column 9).

    The energy is binned by J_perp and normalized by particle count from
    the collisional run (same normalization as ISHW 2022 Matlab code).

    Parameters
    ----------
    data_coll : LossData
        Loss statistics from run WITH collisions.
    data_nocoll : LossData, optional
        Loss statistics from run WITHOUT collisions.
    output_path : str or Path, optional
        If provided, save figure to this path.
    show : bool
        If True, display the figure.
    figsize : tuple
        Figure size in inches.
    title : str, optional
        Plot title.
    xlim : float, optional
        X-axis limit for lost energy fraction.
    nperp : int
        Number of J_perp bins (default 100, same as ISHW 2022).
    smooth_sigma : float
        Gaussian smoothing sigma for curves (default 2.0).

    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    if not HAS_PLOTTING:
        raise ImportError("Plotting requires matplotlib and scipy")

    from scipy.ndimage import gaussian_filter1d

    fig, ax = plt.subplots(figsize=figsize)

    def compute_binned_energy(data: LossData, pmax: float, hp: float):
        """Bin particles and energy by J_perp."""
        part_distr = np.zeros(nperp)
        energ_distr = np.zeros(nperp)

        for i in range(data.n_particles):
            k = int(np.ceil(data.perp_invariant[i] / hp)) - 1
            k = min(nperp - 1, max(0, k))
            part_distr[k] += 1

            if data.lost_mask[i]:
                energ_distr[k] += data.final_p[i] ** 2

        return part_distr, energ_distr

    pmax = np.max(data_coll.perp_invariant)
    if data_nocoll is not None:
        pmax = max(pmax, np.max(data_nocoll.perp_invariant))
    hp = pmax / nperp

    part_c, energ_c = compute_binned_energy(data_coll, pmax, hp)

    if data_nocoll is not None:
        _, energ_n = compute_binned_energy(data_nocoll, pmax, hp)
    else:
        energ_n = None

    part_c_safe = np.where(part_c > 0, part_c, 1)
    jperp = np.arange(1, nperp + 1) / nperp

    enl_coll = gaussian_filter1d(energ_c / part_c_safe, smooth_sigma)

    max_f = np.max(enl_coll)

    if energ_n is not None:
        enl_nocoll = gaussian_filter1d(energ_n / part_c_safe, smooth_sigma)
        ax.plot(enl_nocoll, jperp, "b-", lw=2, label="No coll")
        max_f = max(max_f, np.max(enl_nocoll))

    ax.plot(enl_coll, jperp, "r-", lw=2, label="With coll")

    ax.set_ylim([0, 1])
    if xlim is not None:
        ax.set_xlim([0, xlim])
    else:
        ax.set_xlim([0, max_f * 1.1])

    ax.set_xlabel("lost energy fraction")
    ax.set_ylabel(r"$J_\perp$")
    ax.legend()

    if title:
        ax.set_title(title)

    plt.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_kde_loss_density(
    data: LossData,
    output_path: Optional[str | Path] = None,
    show: bool = True,
    figsize: Tuple[float, float] = (5, 4),
) -> Figure:
    """
    Plot KDE density of loss time vs trapping parameter (Brazilian flag plot).

    This is the signature visualization from the ISHW/Accelerated Methods papers,
    showing where in parameter space particles are lost.

    Parameters
    ----------
    data : LossData
        Loss statistics data.
    output_path : str or Path, optional
        If provided, save figure to this path.
    show : bool
        If True, display the figure.
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    if not HAS_PLOTTING:
        raise ImportError("Plotting requires matplotlib and scipy")

    fig, ax = plt.subplots(figsize=figsize)

    lost = data.lost_mask
    tlost = np.abs(data.loss_times[lost])
    trap_par = data.trap_parameter[lost]

    t_threshold = 0.99 * data.trace_time
    kde_mask = tlost < t_threshold

    if np.sum(kde_mask) > 10:
        log_tlost = np.log10(tlost[kde_mask])
        trap_kde = trap_par[kde_mask]

        try:
            kde = gaussian_kde([log_tlost, trap_kde])
            X, Y = np.mgrid[-5:0:100j, -2:1:150j]
            positions = np.vstack([X.ravel(), Y.ravel()])
            Z = np.reshape(kde(positions).T, X.shape)

            ax.imshow(
                np.rot90(Z),
                cmap=plt.cm.Blues,
                extent=[-5, 0, -2, 1],
                aspect="auto",
            )
        except np.linalg.LinAlgError:
            pass

    if np.sum(lost) > 0:
        ax.plot(np.log10(tlost), trap_par, "k,", markersize=0.5, alpha=0.5)

    ax.axhline(y=0, color="k", linestyle="--", linewidth=0.8)
    ax.set_ylim([-2.1, 1.0])
    ax.set_xlim([-5.0, 0.1])
    ax.set_ylabel("Trapping parameter $\\eta$")
    ax.set_xlabel(r"Loss time $t$ / s")
    ax.set_xticks([-5, -4, -3, -2, -1, 0])
    ax.set_xticklabels(
        [r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$1$"]
    )

    # Twin axis for confined fraction
    ax2 = ax.twinx()
    t_conf = data.time_grid
    valid_t = t_conf > 0
    ax2.plot(np.log10(t_conf[valid_t]), data.confined_fraction[valid_t], color="tab:red")
    ax2.set_ylim([0.5, 1.25])
    ax2.set_ylabel("Confined fraction", color="tab:red")
    ax2.tick_params("y", colors="tab:red")

    plt.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig
