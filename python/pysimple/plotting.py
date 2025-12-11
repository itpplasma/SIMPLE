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

# Physics constants for alpha particle slowing down in DT plasma
ALPHA_BIRTH_ENERGY_EV = 3.5e6  # eV - DT fusion alpha particle birth energy
THERMAL_ENERGY_EV = 1.5e4  # eV - typical thermal energy at 10 keV
THERMAL_ENERGY_FRACTION = THERMAL_ENERGY_EV / ALPHA_BIRTH_ENERGY_EV


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
    slowing_down_curve: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute energy confined fraction over time.

    The energy confined fraction accounts for the fact that particles
    lost early carry away more energy than those lost late (after
    slowing down via collisions).

    Parameters
    ----------
    data : LossData
        Loss statistics data from load_loss_data().
    time_grid : np.ndarray, optional
        Time points to evaluate. Defaults to data.time_grid.
    slowing_down_curve : tuple of (times, energy), optional
        If provided, use theoretical energy from slowing-down curve lookup.
        Otherwise use actual final_p^2 from simulation.

    Returns
    -------
    times : np.ndarray
        Time grid [s].
    energy_fraction : np.ndarray
        Energy confined fraction at each time (0 to 1).
    """
    if time_grid is None:
        time_grid = data.time_grid

    e_total = float(data.n_particles)

    if slowing_down_curve is not None:
        times_sd, energy_sd = slowing_down_curve
        energy_at_loss = np.where(
            data.lost_mask,
            np.interp(data.loss_times, times_sd, energy_sd, left=energy_sd[0], right=energy_sd[-1]),
            0.0,
        )
    else:
        energy_at_loss = data.final_p**2

    energy_fraction = np.zeros_like(time_grid)
    for i, t in enumerate(time_grid):
        lost_before_t = (data.loss_times > 0) & (data.loss_times < t)
        e_lost = np.sum(energy_at_loss[lost_before_t])
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


def _plot_confined_fraction_panel(ax, data: LossData) -> None:
    """Plot confined fraction vs time panel (particles and energy)."""
    t_conf = data.time_grid
    conf_frac = data.confined_fraction
    t_energy, energy_frac = compute_energy_confined_fraction(data)

    valid_t = t_conf > 0
    ax.semilogx(t_conf[valid_t], conf_frac[valid_t], "r-", linewidth=1.5, label="Particles")
    ax.semilogx(t_energy[valid_t], energy_frac[valid_t], "b-", linewidth=1.5, label="Energy")

    ax.set_xlabel(r"Time $t$ / s")
    ax.set_ylabel("Confined fraction")
    ax.set_ylim([0, 1.05])
    ax.legend(loc="lower left")
    ax.grid(True, alpha=0.3)
    ax.set_title("(a) Confined fraction vs time")


def _plot_kde_density_panel(ax, data: LossData) -> None:
    """Plot KDE density of loss time vs trapping parameter panel."""
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
            ax.imshow(np.rot90(Z), cmap=plt.cm.Blues, extent=[-5, 0, -2, 1], aspect="auto")
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
    ax.set_xticklabels([r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$1$"])

    ax_twin = ax.twinx()
    t_conf = data.time_grid
    valid_t = t_conf > 0
    ax_twin.plot(np.log10(t_conf[valid_t]), data.confined_fraction[valid_t], color="tab:red", lw=1.5)
    ax_twin.set_ylim([0.5, 1.25])
    ax_twin.set_ylabel("Confined fraction", color="tab:red")
    ax_twin.tick_params("y", colors="tab:red")
    ax.set_title("(b) Loss time vs trapping parameter (KDE)")


def _plot_energy_loss_panel(ax, data: LossData) -> None:
    """Plot energy loss distribution vs J_perp panel."""
    jperp_bins, particle_count, energy_lost = compute_energy_loss_distribution(data)
    total = data.n_particles if data.n_particles > 0 else 1
    energy_frac = energy_lost / total
    particle_frac = particle_count / total

    ax.bar(jperp_bins, energy_frac, width=1.0 / len(jperp_bins), alpha=0.7, label="Energy", color="tab:blue")
    ax.bar(jperp_bins, particle_frac, width=1.0 / len(jperp_bins), alpha=0.4, label="Particles", color="tab:red")
    ax.set_xlabel(r"$J_\perp$ (normalized)")
    ax.set_ylabel("Lost fraction per bin")
    ax.legend()
    ax.set_title("(c) Loss distribution vs $J_\\perp$")


def _plot_starting_positions_panel(ax, data: LossData) -> None:
    """Plot starting positions colored by loss time panel."""
    tlost_all = np.maximum(np.abs(data.loss_times), 1e-10)
    scatter = ax.scatter(data.start_theta, data.start_pitch, c=np.log10(tlost_all), s=1, cmap="viridis", vmin=-5, vmax=0)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(r"$\log_{10}(t_{\mathrm{loss}})$")
    ax.set_xlabel(r"$\theta$ (poloidal angle)")
    ax.set_ylabel(r"$\lambda = v_\parallel / v$ (pitch)")
    ax.set_xlim([0, 2 * np.pi])
    ax.set_ylim([-1, 1])
    ax.set_title("(d) Starting positions colored by loss time")


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
    - (a) Confined fraction vs time (particles and energy)
    - (b) KDE density: loss time vs trapping parameter
    - (c) Energy loss distribution vs J_perp
    - (d) Starting positions colored by loss time

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
    """
    if not HAS_PLOTTING:
        raise ImportError("Plotting requires matplotlib and scipy.")

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    _plot_confined_fraction_panel(axes[0, 0], data)
    _plot_kde_density_panel(axes[0, 1], data)
    _plot_energy_loss_panel(axes[1, 0], data)
    _plot_starting_positions_panel(axes[1, 1], data)

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
    slowing_down_curve: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (7, 5),
) -> Figure:
    """
    Plot confined fraction vs time for particles and energy.

    Supports both actual energy (from final_p^2) and theoretical energy
    (from slowing-down curve lookup).

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
    slowing_down_curve : tuple of (times, energy), optional
        If provided along with include_energy=True, plots both actual
        and theoretical energy confined fractions.
    title : str, optional
        Plot title.
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    if not HAS_PLOTTING:
        raise ImportError("Plotting requires matplotlib")

    fig, ax = plt.subplots(figsize=figsize)

    t_conf = data.time_grid
    valid_t = t_conf > 0

    ax.semilogx(
        t_conf[valid_t],
        data.confined_fraction[valid_t],
        "k-",
        linewidth=2,
        label="Particles",
    )

    if include_energy:
        t_energy, energy_frac = compute_energy_confined_fraction(data)
        ax.semilogx(
            t_energy[valid_t],
            energy_frac[valid_t],
            "b-",
            linewidth=2,
            label="Energy (actual)",
        )

        if slowing_down_curve is not None:
            _, energy_frac_theo = compute_energy_confined_fraction(
                data, slowing_down_curve=slowing_down_curve
            )
            ax.semilogx(
                t_energy[valid_t],
                energy_frac_theo[valid_t],
                "b--",
                linewidth=1.5,
                label="Energy (theoretical)",
            )

    ax.set_xlabel(r"Time $t$ / s")
    ax.set_ylabel("Confined fraction")
    ax.set_ylim([0, 1.05])
    ax.legend(loc="lower left")
    ax.grid(True, alpha=0.3)

    if title:
        ax.set_title(title)

    plt.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def get_slowing_down_curve_path() -> Path:
    """
    Get path to bundled slowing-down energy curve data file.

    Returns the path to energyslow_aver.dat included with pysimple.
    This file contains the theoretical alpha particle slowing-down
    curve for 3.5 MeV alphas in a DT plasma.

    Returns
    -------
    Path
        Path to the bundled energyslow_aver.dat file.
    """
    import importlib.resources

    try:
        # Python 3.9+
        return importlib.resources.files("pysimple").joinpath("data/energyslow_aver.dat")
    except AttributeError:
        # Fallback for older Python
        import pkg_resources

        return Path(
            pkg_resources.resource_filename("pysimple", "data/energyslow_aver.dat")
        )


def load_slowing_down_curve(
    filepath: Optional[str | Path] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load theoretical slowing-down energy curve.

    The file contains time vs normalized energy (E/E0) for alpha particles
    slowing down in a DT plasma over one slowing-down time.

    Parameters
    ----------
    filepath : str or Path, optional
        Path to energyslow_aver.dat file. If None, uses the bundled
        data file included with pysimple.

    Returns
    -------
    times : np.ndarray
        Time grid [s], normalized to slowing-down time (0 to 1).
    energy : np.ndarray
        Normalized energy at each time (1 at t=0, decreasing to thermal).
    """
    if filepath is None:
        filepath = get_slowing_down_curve_path()
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]


def _compute_bin_indices(perp_invariant: np.ndarray, hp: float, nperp: int) -> np.ndarray:
    """Compute J_perp bin indices for all particles (vectorized)."""
    k = np.ceil(perp_invariant / hp).astype(int) - 1
    return np.clip(k, 0, nperp - 1)


def _bin_energy_actual(
    data: LossData, hp: float, nperp: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Bin particles and energy by J_perp using actual final_p^2 (vectorized)."""
    k = _compute_bin_indices(data.perp_invariant, hp, nperp)
    part_distr = np.bincount(k, minlength=nperp).astype(float)
    energy_weights = np.where(data.lost_mask, data.final_p**2, 0.0)
    energ_distr = np.bincount(k, weights=energy_weights, minlength=nperp)
    return part_distr, energ_distr


def _bin_energy_theoretical(
    data: LossData,
    hp: float,
    nperp: int,
    times_sd: np.ndarray,
    energy_sd: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Bin particles and energy by J_perp using theoretical slowing-down (vectorized)."""
    k = _compute_bin_indices(data.perp_invariant, hp, nperp)
    part_distr = np.bincount(k, minlength=nperp).astype(float)

    energy_adj = energy_sd - THERMAL_ENERGY_FRACTION
    energy_at_loss = np.interp(data.loss_times, times_sd, energy_adj, left=energy_adj[0], right=energy_adj[-1])
    energy_weights = np.where(data.lost_mask, energy_at_loss, 0.0)
    energ_distr = np.bincount(k, weights=energy_weights, minlength=nperp)
    return part_distr, energ_distr


def plot_energy_loss_vs_jperp(
    data_coll: LossData,
    data_nocoll: Optional[LossData] = None,
    output_path: Optional[str | Path] = None,
    show: bool = True,
    figsize: Tuple[float, float] = (7, 6),
    title: Optional[str] = None,
    xlim: Optional[float] = None,
    nperp: int = 100,
    slowing_down_curve: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> Figure:
    """
    Plot lost energy fraction vs J_perp in ISHW 2022 style.

    Supports up to four curves:
    1. COLL (actual): Energy from final_p^2 (actual velocity at loss)
    2. NOCOLL (actual): Energy from final_p^2 (should be ~1 without collisions)
    3. COLL (theoretical): Energy from slowing-down curve lookup by loss time
    4. NOCOLL (theoretical): Energy from slowing-down curve lookup by loss time

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
    slowing_down_curve : tuple of (times, energy), optional
        Theoretical slowing-down curve from load_slowing_down_curve().
        If provided, adds two additional curves using theoretical energy
        lookup instead of actual final_p^2.

    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    if not HAS_PLOTTING:
        raise ImportError("Plotting requires matplotlib and scipy")

    fig, ax = plt.subplots(figsize=figsize)

    pmax = np.max(data_coll.perp_invariant)
    if data_nocoll is not None:
        pmax = max(pmax, np.max(data_nocoll.perp_invariant))
    hp = pmax / nperp

    part_c, energ_c = _bin_energy_actual(data_coll, hp, nperp)
    part_n, energ_n = _bin_energy_actual(data_nocoll, hp, nperp) if data_nocoll else (None, None)

    energ_c_theo, energ_n_theo = None, None
    if slowing_down_curve is not None:
        times_sd, energy_sd = slowing_down_curve
        _, energ_c_theo = _bin_energy_theoretical(data_coll, hp, nperp, times_sd, energy_sd)
        if data_nocoll is not None:
            _, energ_n_theo = _bin_energy_theoretical(data_nocoll, hp, nperp, times_sd, energy_sd)

    part_c_safe = np.where(part_c > 0, part_c, 1)
    jperp = np.arange(1, nperp + 1) / nperp

    curves = [
        (energ_c / part_c_safe, "r-", 2, "COLL (actual)"),
        (energ_n / part_c_safe if energ_n is not None else None, "b-", 2, "NOCOLL (actual)"),
        (energ_c_theo / part_c_safe if energ_c_theo is not None else None, "r--", 1.5, "COLL (theoretical)"),
        (energ_n_theo / part_c_safe if energ_n_theo is not None else None, "b--", 1.5, "NOCOLL (theoretical)"),
    ]

    max_f = 0.0
    for y_data, style, lw, label in curves:
        if y_data is not None:
            ax.plot(y_data, jperp, style, lw=lw, label=label)
            max_f = max(max_f, np.max(y_data))

    ax.set_ylim([0, 1])
    ax.set_xlim([0, xlim if xlim is not None else max_f * 1.1])
    ax.set_xlabel("lost energy fraction")
    ax.set_ylabel(r"$J_\perp$")
    ax.legend(loc="lower right")
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
