"""
Engineering tools for wall heat load analysis and port optimization.

This module provides:
- WallHeatMap: Bin particle wall hits to heat flux (MW/m^2)
- PortOptimizer: Optimize port placement to minimize heat exposure
- Visualization: 3D heat flux plots via PyVista

Physics Background
------------------
For D-T fusion, each reaction produces 17.6 MeV total energy:
- 14.1 MeV -> neutron (escapes)
- 3.5 MeV -> alpha particle (heats plasma or hits wall)

The total alpha power for a reactor is:
    P_alpha = P_fusion / 5

For example, a 3 GW fusion plant has ~600 MW of alpha power.

Heat flux calculation with energy weighting
-------------------------------------------
The code accounts for collisional slowing down of alpha particles.
Each particle's final energy is:

    E_i = p_i^2 * E_birth

where p_i = v/v0 is the normalized velocity at wall impact (zend[3]).

Without collisions: p = 1, E_final = E_birth
With collisions: p < 1, E_final < E_birth (particles slowed down)

Heat flux per bin:

    q_bin = sum(p_i^2 for i in bin) * (P_alpha / N_total) / A_bin

Lost energy fraction:

    f_energy = sum(p_i^2 for lost particles) / N_total

This is more accurate than particle-based loss fraction because:
- Prompt losses (early escape, p ~ 1) deposit more energy
- Slow losses (after slowing down, p << 1) deposit less energy

Typical values (from literature):
- Infinity Two: ~2.5 MW/m^2 peak
- ARIES-CS: 5-18 MW/m^2 peak
- W7-X experiments: 0.1-1 MW/m^2 (fast ion loads)

References:
- Lazerson et al., Plasma Phys. Control. Fusion 63 (2021) 125033
- Albert et al., J. Plasma Phys. (2024) - Infinity Two study
- Ku & Boozer, Fusion Sci. Technol. 54 (2008) 673 - ARIES-CS

Example:
    from pysimple.engineering import WallHeatMap, PortOptimizer

    # For a 3 GW fusion reactor (600 MW alpha power)
    heat_map = WallHeatMap.from_netcdf(
        "results.nc",
        total_alpha_power_MW=600.0,  # Required for physical MW/m^2
        major_radius=1000.0,  # cm
        minor_radius=100.0,   # cm
    )
    print(f"Loss fraction: {heat_map.loss_fraction:.1%}")
    print(f"Peak flux: {heat_map.peak_flux:.2f} MW/m^2")

    opt = PortOptimizer(heat_map)
    opt.add_port("NBI", theta_width=0.3, zeta_width=0.2)
    result = opt.solve()
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from pysimple.coil_clearance import (
    CoilClearanceConstraint,
    CoilSegments,
    WallSurface,
    clearance_map_on_heatmap_grid,
)

try:
    import netCDF4 as nc
    HAS_NETCDF4 = True
except ImportError:
    HAS_NETCDF4 = False


def _compute_surface_element_grid(
    x_wall: np.ndarray,
    y_wall: np.ndarray,
    z_wall: np.ndarray,
    dtheta: float,
    dzeta: float,
) -> np.ndarray:
    """
    Compute surface element |dr/dtheta x dr/dzeta| at each grid point.

    Uses central differences for theta (periodic) and zeta interior points.
    Uses one-sided differences at zeta boundaries because Cartesian
    coordinates are not periodic within one field period (they rotate).

    Args:
        x_wall, y_wall, z_wall: Wall coordinates, shape (nzeta, ntheta), in meters
        dtheta: Grid spacing in theta
        dzeta: Grid spacing in zeta

    Returns:
        surface_element: Array of shape (nzeta, ntheta) with |dr/dtheta x dr/dzeta|
    """
    nzeta, ntheta = x_wall.shape
    surface_element = np.zeros((nzeta, ntheta))

    for iz in range(nzeta):
        for it in range(ntheta):
            # Theta: periodic, use central difference with wrapping
            it_p = (it + 1) % ntheta
            it_m = (it - 1) % ntheta
            dr_dtheta = np.array([
                (x_wall[iz, it_p] - x_wall[iz, it_m]) / (2 * dtheta),
                (y_wall[iz, it_p] - y_wall[iz, it_m]) / (2 * dtheta),
                (z_wall[iz, it_p] - z_wall[iz, it_m]) / (2 * dtheta),
            ])

            # Zeta: one-sided differences at boundaries
            # (Cartesian coords rotate by 2*pi/nfp, not periodic in one period)
            if iz == 0:
                dr_dzeta = np.array([
                    (x_wall[1, it] - x_wall[0, it]) / dzeta,
                    (y_wall[1, it] - y_wall[0, it]) / dzeta,
                    (z_wall[1, it] - z_wall[0, it]) / dzeta,
                ])
            elif iz == nzeta - 1:
                dr_dzeta = np.array([
                    (x_wall[iz, it] - x_wall[iz - 1, it]) / dzeta,
                    (y_wall[iz, it] - y_wall[iz - 1, it]) / dzeta,
                    (z_wall[iz, it] - z_wall[iz - 1, it]) / dzeta,
                ])
            else:
                dr_dzeta = np.array([
                    (x_wall[iz + 1, it] - x_wall[iz - 1, it]) / (2 * dzeta),
                    (y_wall[iz + 1, it] - y_wall[iz - 1, it]) / (2 * dzeta),
                    (z_wall[iz + 1, it] - z_wall[iz - 1, it]) / (2 * dzeta),
                ])

            cross = np.cross(dr_dtheta, dr_dzeta)
            surface_element[iz, it] = np.sqrt(np.sum(cross**2))

    return surface_element


def compute_wall_area_from_chartmap(chartmap_file: str | Path) -> float:
    """
    Compute actual wall surface area from chartmap NetCDF file.

    The chartmap contains Cartesian coordinates (x, y, z) on a grid of
    (rho, theta, zeta). The wall surface is at rho=1 (last index).

    Surface area is computed as:
        A_per_period = integral |dr/dtheta x dr/dzeta| dtheta dzeta
        A_total = n_fp * A_per_period

    Args:
        chartmap_file: Path to chartmap NetCDF file

    Returns:
        Total wall surface area in m^2

    Note:
        This gives the actual 3D shaped stellarator wall area, which can be
        1.5-2x larger than the simple torus approximation 4*pi^2*R*a.
    """
    if not HAS_NETCDF4:
        raise ImportError("netCDF4 required: pip install netCDF4")

    with nc.Dataset(chartmap_file, 'r') as ds:
        theta = ds.variables['theta'][:]
        zeta = ds.variables['zeta'][:]
        nfp = int(ds.variables['num_field_periods'][...])

        x_wall = ds.variables['x'][:, :, -1] * 0.01  # cm -> m
        y_wall = ds.variables['y'][:, :, -1] * 0.01
        z_wall = ds.variables['z'][:, :, -1] * 0.01

    nzeta, ntheta = x_wall.shape
    dtheta = theta[1] - theta[0] if len(theta) > 1 else 2 * np.pi
    dzeta = zeta[1] - zeta[0] if len(zeta) > 1 else 2 * np.pi / nfp

    surface_element = _compute_surface_element_grid(
        x_wall, y_wall, z_wall, dtheta, dzeta
    )

    return float(np.sum(surface_element) * dtheta * dzeta * nfp)


def compute_bin_areas_from_chartmap(
    chartmap_file: str | Path,
    theta_edges: np.ndarray,
    zeta_edges: np.ndarray,
) -> np.ndarray:
    """
    Compute actual area of each (theta, zeta) bin from chartmap geometry.

    This properly accounts for the 3D stellarator shape rather than using
    a uniform area approximation.

    Coordinate handling:
        - Chartmap theta: typically [0, 2π], mapped to heat map [-π, π]
        - Chartmap zeta: one field period [0, 2π/nfp], extended by symmetry
        - Heat map zeta: full torus [0, 2π], each bin maps to one period

    Args:
        chartmap_file: Path to chartmap NetCDF file
        theta_edges: Bin edges for poloidal angle (n_theta + 1), typically [-π, π]
        zeta_edges: Bin edges for toroidal angle (n_zeta + 1), typically [0, 2π]

    Returns:
        bin_areas: Array of shape (n_theta, n_zeta) with area in m^2 per bin
    """
    if not HAS_NETCDF4:
        raise ImportError("netCDF4 required: pip install netCDF4")

    with nc.Dataset(chartmap_file, 'r') as ds:
        theta_cm = ds.variables['theta'][:]
        zeta_cm = ds.variables['zeta'][:]
        nfp = int(ds.variables['num_field_periods'][...])

        x_wall = ds.variables['x'][:, :, -1] * 0.01  # cm -> m
        y_wall = ds.variables['y'][:, :, -1] * 0.01
        z_wall = ds.variables['z'][:, :, -1] * 0.01

    nzeta_cm, ntheta_cm = x_wall.shape
    dtheta_cm = theta_cm[1] - theta_cm[0] if len(theta_cm) > 1 else 2 * np.pi
    dzeta_cm = zeta_cm[1] - zeta_cm[0] if len(zeta_cm) > 1 else 2 * np.pi / nfp

    # Compute surface element using shared helper (consistent boundary handling)
    surface_element = _compute_surface_element_grid(
        x_wall, y_wall, z_wall, dtheta_cm, dzeta_cm
    )

    # Convert chartmap theta [0, 2π] to heat map convention [-π, π]
    # theta in [0, π] stays same, theta in (π, 2π] -> (-π, 0]
    theta_hm = np.where(theta_cm > np.pi, theta_cm - 2 * np.pi, theta_cm)

    # Zeta period for the chartmap
    zeta_period = 2 * np.pi / nfp

    # Integrate surface element over each output bin
    n_theta = len(theta_edges) - 1
    n_zeta = len(zeta_edges) - 1
    bin_areas = np.zeros((n_theta, n_zeta))

    for i_theta in range(n_theta):
        th_min, th_max = theta_edges[i_theta], theta_edges[i_theta + 1]
        for i_zeta in range(n_zeta):
            ze_min, ze_max = zeta_edges[i_zeta], zeta_edges[i_zeta + 1]

            # Map heat map zeta bin to chartmap period [0, zeta_period]
            # Due to stellarator symmetry, all field periods have same geometry
            ze_min_cm = ze_min % zeta_period
            ze_max_cm = ze_max % zeta_period

            # Find chartmap grid points in this bin and sum their contributions
            area_sum = 0.0
            for iz in range(nzeta_cm):
                zeta_val = zeta_cm[iz]
                # Handle wrapped bins (ze_min_cm > ze_max_cm after modulo)
                if ze_min_cm <= ze_max_cm:
                    in_zeta_bin = ze_min_cm <= zeta_val < ze_max_cm
                else:
                    # Bin wraps around: includes [ze_min_cm, period) and [0, ze_max_cm)
                    in_zeta_bin = zeta_val >= ze_min_cm or zeta_val < ze_max_cm

                if not in_zeta_bin:
                    continue
                for it in range(ntheta_cm):
                    theta_val = theta_hm[it]
                    if not (th_min <= theta_val < th_max):
                        continue
                    area_sum += surface_element[iz, it] * dtheta_cm * dzeta_cm

            bin_areas[i_theta, i_zeta] = area_sum

    return bin_areas


@dataclass
class WallHeatMap:
    """
    Wall heat flux distribution from SIMPLE particle tracing results.

    Bins lost particles by their wall hit location (theta, zeta) and computes
    heat flux in MW/m^2. Accounts for collisional slowing down by using
    energy-weighted binning based on final velocity p = v/v0.

    The heat flux calculation is:
        q_bin = sum(p_i^2 for i in bin) * (P_alpha / N_total) / A_bin

    where:
        - p_i = v/v0: normalized final velocity (from zend[3])
        - P_alpha: total alpha power (MW), user-provided
        - N_total: total particles traced
        - A_bin: area of this bin (m^2), from chartmap geometry or uniform approx

    Coordinate Conventions (VMEC):
        - theta: Poloidal angle in radians, range [-pi, pi]
            * theta = 0: OUTBOARD midplane (maximum R, Z=0)
            * theta = +/- pi: INBOARD midplane (minimum R, Z=0)
            * theta > 0: Upper half (Z > 0)
            * theta < 0: Lower half (Z < 0)
        - zeta: Toroidal angle in radians, range [0, 2*pi]
            * zeta increases in the direction of the toroidal current
            * For nfp field periods, geometry repeats every 2*pi/nfp

    For port placement:
        - Outboard side: |theta| < pi/2 (e.g., NBI ports)
        - Inboard side: |theta| > pi/2

    Attributes:
        theta_edges: Bin edges for poloidal angle (rad)
        zeta_edges: Bin edges for toroidal angle (rad)
        flux_grid: Heat flux in MW/m^2, shape (n_theta, n_zeta)
        energy_grid: Energy weight sum per bin (sum of p^2)
        hit_count: Number of particles per bin
        lost_power: Power lost to wall (MW) = energy_loss_fraction * total_alpha_power
        total_alpha_power: Total alpha power from fusion (MW), user input
        n_lost: Number of lost particles
        n_total: Total number of particles traced
        trace_time: Tracing time in seconds
        wall_area: Total wall area in m^2
        energy_loss_fraction: Fraction of energy lost = sum(p^2) / N_total
        bin_areas: Per-bin areas in m^2, shape (n_theta, n_zeta)
        flux_error: Monte Carlo error estimate (MW/m^2), = flux / sqrt(N_bin)
    """

    theta_edges: np.ndarray
    zeta_edges: np.ndarray
    flux_grid: np.ndarray
    energy_grid: np.ndarray
    hit_count: np.ndarray
    lost_power: float
    total_alpha_power: float
    energy_loss_fraction: float
    n_lost: int
    n_total: int
    trace_time: float
    wall_area: float = 0.0
    major_radius: float = 0.0
    minor_radius: float = 0.0
    bin_areas: Optional[np.ndarray] = field(default=None, repr=False)
    flux_error: Optional[np.ndarray] = field(default=None, repr=False)
    _wall_positions: Optional[np.ndarray] = field(default=None, repr=False)
    _loss_times: Optional[np.ndarray] = field(default=None, repr=False)
    _final_p: Optional[np.ndarray] = field(default=None, repr=False)

    @classmethod
    def from_netcdf(
        cls,
        results_file: str | Path,
        total_alpha_power_MW: float,
        n_theta: int = 64,
        n_zeta: int = 128,
        major_radius: float = 1000.0,
        minor_radius: float = 100.0,
        chartmap_file: Optional[str | Path] = None,
    ) -> "WallHeatMap":
        """
        Load results from NetCDF and compute heat flux distribution.

        The heat flux is computed using the Monte Carlo approach:
            q_bin = sum(p_i^2 for i in bin) * (P_alpha / N_total) / A_bin

        This requires the user to specify total_alpha_power_MW based on the
        plasma scenario. For D-T fusion: P_alpha = P_fusion / 5.

        Args:
            results_file: Path to results.nc from SIMPLE
            total_alpha_power_MW: Total alpha power from fusion reactions (MW).
                For D-T: P_alpha = P_fusion / 5. E.g., 3 GW fusion -> 600 MW.
            n_theta: Number of poloidal bins
            n_zeta: Number of toroidal bins
            major_radius: Major radius in cm (for torus wall area approximation)
            minor_radius: Minor radius in cm (for torus wall area approximation)
            chartmap_file: Optional path to chartmap NetCDF for actual wall geometry.
                If provided, uses real 3D surface area from metric tensor instead
                of simple torus approximation.

        Returns:
            WallHeatMap with physically-scaled heat flux in MW/m^2

        Example:
            # With actual wall geometry from chartmap
            heat_map = WallHeatMap.from_netcdf(
                "results.nc",
                total_alpha_power_MW=600.0,
                chartmap_file="wout.chartmap.nc",
            )

        Note:
            Guiding center limitation: The reported wall positions are guiding
            center positions, not actual wall impact points. For 3.5 MeV alphas,
            the gyroradius is ~5-10 cm, which smears the heat load pattern.
        """
        if not HAS_NETCDF4:
            raise ImportError("netCDF4 required: pip install netCDF4")

        if total_alpha_power_MW <= 0:
            raise ValueError(
                f"total_alpha_power_MW must be positive, got {total_alpha_power_MW}"
            )

        results_file = Path(results_file)
        if not results_file.exists():
            raise FileNotFoundError(f"Results file not found: {results_file}")

        with nc.Dataset(results_file, 'r') as ds:
            zend = ds.variables['zend'][:]
            xend_cart = ds.variables['xend_cart'][:]
            class_lost = ds.variables['class_lost'][:].astype(bool)
            times_lost = ds.variables['times_lost'][:]
            trace_time = float(ds.trace_time)
            n_total = int(ds.ntestpart)

        lost_mask = class_lost
        n_lost = int(np.sum(lost_mask))

        theta_edges = np.linspace(-np.pi, np.pi, n_theta + 1)
        zeta_edges = np.linspace(0, 2 * np.pi, n_zeta + 1)

        # Compute wall area: use chartmap if provided, else torus approximation
        if chartmap_file is not None:
            wall_area_m2 = compute_wall_area_from_chartmap(chartmap_file)
        else:
            # Simple torus: A = 4 * pi^2 * R0 * a, convert cm to m
            wall_area_m2 = 4 * np.pi**2 * major_radius * minor_radius * 1e-4

        if n_lost == 0:
            return cls(
                theta_edges=theta_edges,
                zeta_edges=zeta_edges,
                flux_grid=np.zeros((n_theta, n_zeta)),
                energy_grid=np.zeros((n_theta, n_zeta)),
                hit_count=np.zeros((n_theta, n_zeta), dtype=int),
                lost_power=0.0,
                total_alpha_power=total_alpha_power_MW,
                energy_loss_fraction=0.0,
                n_lost=0,
                n_total=n_total,
                trace_time=trace_time,
                wall_area=wall_area_m2,
                major_radius=major_radius,
                minor_radius=minor_radius,
            )

        theta_lost = zend[1, lost_mask]
        zeta_lost = zend[2, lost_mask]
        # zend[3] = p = v/v0 = normalized velocity at wall impact
        final_p = zend[3, lost_mask]
        wall_positions = xend_cart[:, lost_mask]
        loss_times = times_lost[lost_mask]

        # Map angles to canonical ranges for histogram binning
        # SIMPLE outputs theta in VMEC convention [0, 2pi], map to [-pi, pi]
        # SIMPLE outputs zeta as cumulative angle (can be large), map to [0, 2pi]
        theta_lost = ((theta_lost + np.pi) % (2 * np.pi)) - np.pi
        zeta_lost = zeta_lost % (2 * np.pi)

        # Energy weight = p^2 (normalized to birth energy)
        energy_weight = final_p**2

        # Bin by particle count
        hit_count, _, _ = np.histogram2d(
            theta_lost, zeta_lost,
            bins=[theta_edges, zeta_edges]
        )
        hit_count = hit_count.astype(int)

        # Energy-weighted histogram: sum of p^2 per bin
        energy_grid, _, _ = np.histogram2d(
            theta_lost, zeta_lost,
            bins=[theta_edges, zeta_edges],
            weights=energy_weight
        )

        # Compute bin areas: use 3D geometry if chartmap provided, else uniform
        if chartmap_file is not None:
            bin_areas = compute_bin_areas_from_chartmap(
                chartmap_file, theta_edges, zeta_edges
            )
            # Warn if bin resolution exceeds chartmap resolution
            bin_area_sum = np.sum(bin_areas)
            if bin_area_sum < 0.9 * wall_area_m2:
                import warnings
                warnings.warn(
                    f"Bin areas sum ({bin_area_sum:.1f} m^2) is less than 90% of "
                    f"wall area ({wall_area_m2:.1f} m^2). Heat map resolution may "
                    f"exceed chartmap resolution, causing empty bins. Consider "
                    f"reducing n_theta/n_zeta or using a finer chartmap.",
                    UserWarning,
                    stacklevel=2,
                )
        else:
            # Uniform bin area approximation (valid for simple torus)
            uniform_area = wall_area_m2 / (n_theta * n_zeta)
            bin_areas = np.full((n_theta, n_zeta), uniform_area)

        # Energy-weighted heat flux calculation:
        # q_bin = sum(p_i^2 for i in bin) * (P_alpha / n_total) / A_bin
        power_density = total_alpha_power_MW / n_total  # MW per MC particle
        # Avoid division by zero for empty bins
        with np.errstate(divide='ignore', invalid='ignore'):
            flux_grid = np.where(
                bin_areas > 0,
                energy_grid * power_density / bin_areas,
                0.0
            )

        # Monte Carlo error estimate: sigma_q = q / sqrt(N_bin)
        # For bins with few particles, error is large
        with np.errstate(divide='ignore', invalid='ignore'):
            flux_error = np.where(
                hit_count > 0,
                flux_grid / np.sqrt(hit_count),
                0.0
            )

        # Energy loss fraction = sum(p^2 for all lost) / n_total
        total_energy_lost = np.sum(energy_weight)
        energy_loss_fraction = total_energy_lost / n_total

        # Lost power based on energy, not particle count
        lost_power = energy_loss_fraction * total_alpha_power_MW

        return cls(
            theta_edges=theta_edges,
            zeta_edges=zeta_edges,
            flux_grid=flux_grid,
            energy_grid=energy_grid,
            hit_count=hit_count,
            lost_power=lost_power,
            total_alpha_power=total_alpha_power_MW,
            energy_loss_fraction=energy_loss_fraction,
            n_lost=n_lost,
            n_total=n_total,
            trace_time=trace_time,
            wall_area=wall_area_m2,
            major_radius=major_radius,
            minor_radius=minor_radius,
            bin_areas=bin_areas,
            flux_error=flux_error,
            _wall_positions=wall_positions,
            _loss_times=loss_times,
            _final_p=final_p,
        )

    @property
    def theta_centers(self) -> np.ndarray:
        """Bin centers for poloidal angle."""
        return 0.5 * (self.theta_edges[:-1] + self.theta_edges[1:])

    @property
    def zeta_centers(self) -> np.ndarray:
        """Bin centers for toroidal angle."""
        return 0.5 * (self.zeta_edges[:-1] + self.zeta_edges[1:])

    @property
    def peak_flux(self) -> float:
        """Maximum heat flux in MW/m^2."""
        return float(np.max(self.flux_grid))

    @property
    def mean_flux(self) -> float:
        """Mean heat flux over non-zero bins in MW/m^2."""
        nonzero = self.flux_grid[self.flux_grid > 0]
        return float(np.mean(nonzero)) if len(nonzero) > 0 else 0.0

    @property
    def loss_fraction(self) -> float:
        """Fraction of particles lost to wall."""
        return self.n_lost / self.n_total if self.n_total > 0 else 0.0

    def flux_at(self, theta: float, zeta: float) -> float:
        """Get heat flux at specific (theta, zeta) location."""
        i_theta = np.searchsorted(self.theta_edges, theta) - 1
        i_zeta = np.searchsorted(self.zeta_edges, zeta) - 1
        i_theta = np.clip(i_theta, 0, len(self.theta_centers) - 1)
        i_zeta = np.clip(i_zeta, 0, len(self.zeta_centers) - 1)
        return float(self.flux_grid[i_theta, i_zeta])

    def integrated_flux_in_region(
        self,
        theta_min: float,
        theta_max: float,
        zeta_min: float,
        zeta_max: float,
    ) -> float:
        """
        Compute total power deposited in a region (MW), energy-weighted.

        Handles zeta wrapping: if zeta_min < 0, splits into two regions
        [2pi + zeta_min, 2pi] and [0, zeta_max].
        """
        theta_mask = (self.theta_centers >= theta_min) & \
                     (self.theta_centers <= theta_max)

        # Handle zeta wrapping for regions near zeta = 0
        if zeta_min < 0:
            # Split into two regions: [2pi + zeta_min, 2pi] and [0, zeta_max]
            zeta_mask = (self.zeta_centers >= 2 * np.pi + zeta_min) | \
                        (self.zeta_centers <= zeta_max)
        elif zeta_max > 2 * np.pi:
            # Split into two regions: [zeta_min, 2pi] and [0, zeta_max - 2pi]
            zeta_mask = (self.zeta_centers >= zeta_min) | \
                        (self.zeta_centers <= zeta_max - 2 * np.pi)
        else:
            zeta_mask = (self.zeta_centers >= zeta_min) & \
                        (self.zeta_centers <= zeta_max)

        # Use energy_grid (sum of p^2) for energy-weighted power
        region_energy = self.energy_grid[np.ix_(theta_mask, zeta_mask)]
        # Power = (sum(p^2) in region / n_total) * P_alpha
        fraction = np.sum(region_energy) / self.n_total if self.n_total > 0 else 0
        return self.total_alpha_power * fraction

    def summary(self) -> str:
        """Return a summary string."""
        lost_pct = 100 * self.lost_power / self.total_alpha_power \
            if self.total_alpha_power > 0 else 0
        return (
            f"WallHeatMap: {self.n_lost}/{self.n_total} particles lost "
            f"({100*self.loss_fraction:.1f}%)\n"
            f"  Alpha power: {self.total_alpha_power:.1f} MW (input)\n"
            f"  Particle loss fraction: {100*self.loss_fraction:.2f}%\n"
            f"  Energy loss fraction: {100*self.energy_loss_fraction:.2f}%\n"
            f"  Lost power: {self.lost_power:.3f} MW ({lost_pct:.1f}%)\n"
            f"  Peak flux: {self.peak_flux:.3f} MW/m^2\n"
            f"  Mean flux: {self.mean_flux:.3f} MW/m^2\n"
            f"  Wall area: {self.wall_area:.1f} m^2"
        )

    def __repr__(self) -> str:
        return self.summary()


@dataclass
class Port:
    """A port or component on the wall surface."""

    name: str
    theta_center: float
    zeta_center: float
    theta_width: float
    zeta_width: float

    @property
    def theta_min(self) -> float:
        return self.theta_center - self.theta_width / 2

    @property
    def theta_max(self) -> float:
        return self.theta_center + self.theta_width / 2

    @property
    def zeta_min(self) -> float:
        return self.zeta_center - self.zeta_width / 2

    @property
    def zeta_max(self) -> float:
        return self.zeta_center + self.zeta_width / 2


@dataclass
class OptimizationResult:
    """Result of port placement optimization."""

    ports: list[Port]
    total_flux_on_ports: float
    max_flux_on_ports: float
    success: bool
    message: str

    def __repr__(self) -> str:
        lines = [f"OptimizationResult(success={self.success})"]
        lines.append(f"  Max flux on ports: {self.max_flux_on_ports:.3f} MW/m^2")
        lines.append(f"  Total flux on ports: {self.total_flux_on_ports:.3f} MW")
        for port in self.ports:
            lines.append(
                f"  {port.name}: theta={port.theta_center:.2f}, "
                f"zeta={port.zeta_center:.2f}"
            )
        return "\n".join(lines)


class PortOptimizer:
    """
    Optimize port placement to minimize heat flux exposure.

    Example:
        opt = PortOptimizer(heat_map)
        opt.add_port("NBI_1", theta_width=0.3, zeta_width=0.2)
        opt.set_max_flux_constraint(0.5)  # MW/m^2
        result = opt.solve()
    """

    def __init__(self, heat_map: WallHeatMap):
        self.heat_map = heat_map
        self._ports: list[dict] = []
        self._exclusion_zones: list[tuple[float, float, float, float]] = []
        self._max_flux_constraint: Optional[float] = None
        self._coil_clearance: Optional[CoilClearanceConstraint] = None

    def add_port(
        self,
        name: str,
        theta_width: float,
        zeta_width: float,
        theta_init: Optional[float] = None,
        zeta_init: Optional[float] = None,
    ) -> "PortOptimizer":
        """
        Add a port to optimize.

        Args:
            name: Port identifier
            theta_width: Poloidal extent in radians
            zeta_width: Toroidal extent in radians
            theta_init: Initial theta position (optional)
            zeta_init: Initial zeta position (optional)
        """
        self._ports.append({
            "name": name,
            "theta_width": theta_width,
            "zeta_width": zeta_width,
            "theta_init": theta_init,
            "zeta_init": zeta_init,
        })
        return self

    def add_exclusion_zone(
        self,
        theta_min: float,
        theta_max: float,
        zeta_min: float = 0.0,
        zeta_max: float = 2 * np.pi,
    ) -> "PortOptimizer":
        """Add a region where ports cannot be placed."""
        self._exclusion_zones.append((theta_min, theta_max, zeta_min, zeta_max))
        return self

    def set_max_flux_constraint(self, max_flux: float) -> "PortOptimizer":
        """Set maximum allowed flux on any port (MW/m^2)."""
        self._max_flux_constraint = max_flux
        return self

    def set_coil_clearance_constraint(
        self,
        coil_file: str | Path,
        *,
        chartmap_file: str | Path,
        min_clearance_m: float,
        use_zero_current_separators: bool = True,
        jump_factor: float = 25.0,
    ) -> "PortOptimizer":
        """
        Enforce a minimum distance from ports to coils (geometric feasibility).

        This constraint requires a chartmap file to map (theta, zeta) port
        locations to physical-space coordinates on the wall.
        """
        if min_clearance_m <= 0.0:
            raise ValueError("min_clearance_m must be positive")

        wall = WallSurface.from_chartmap(chartmap_file)
        coils = CoilSegments.from_file(
            coil_file,
            use_zero_current_separators=use_zero_current_separators,
            jump_factor=jump_factor,
        )
        self._coil_clearance = CoilClearanceConstraint(
            wall=wall,
            coils=coils,
            min_clearance_m=float(min_clearance_m),
        )
        return self

    def _objective(self, x: np.ndarray) -> float:
        """Objective: minimize total flux on all ports, penalize constraints."""
        total = 0.0
        penalty = 0.0
        max_flux_density = 0.0

        for i, port_spec in enumerate(self._ports):
            theta_c = x[2 * i]
            zeta_c = x[2 * i + 1]
            theta_w = port_spec["theta_width"]
            zeta_w = port_spec["zeta_width"]

            theta_min = theta_c - theta_w / 2
            theta_max = theta_c + theta_w / 2
            zeta_min = zeta_c - zeta_w / 2
            zeta_max = zeta_c + zeta_w / 2

            flux = self.heat_map.integrated_flux_in_region(
                theta_min, theta_max, zeta_min, zeta_max
            )
            total += flux

            # Check peak flux density on this port for constraint
            if self._max_flux_constraint is not None:
                port_flux_density = self._get_peak_flux_in_region(
                    theta_min, theta_max, zeta_min, zeta_max
                )
                max_flux_density = max(max_flux_density, port_flux_density)

            # Penalty for port overlapping exclusion zones (check all corners + center)
            if self._port_overlaps_exclusion(theta_c, zeta_c, theta_w, zeta_w):
                penalty += 1e6
            if self._coil_clearance is not None:
                if self._coil_clearance.port_violates(theta_c, zeta_c, theta_w, zeta_w):
                    penalty += 1e6

        # Penalty for exceeding max flux constraint
        if self._max_flux_constraint is not None:
            if max_flux_density > self._max_flux_constraint:
                excess = max_flux_density - self._max_flux_constraint
                penalty += 1e4 * excess  # Proportional penalty

        return total + penalty

    def _get_peak_flux_in_region(
        self, theta_min: float, theta_max: float, zeta_min: float, zeta_max: float
    ) -> float:
        """Get peak flux density (MW/m^2) in a region."""
        theta_mask = (self.heat_map.theta_centers >= theta_min) & \
                     (self.heat_map.theta_centers <= theta_max)
        zeta_mask = (self.heat_map.zeta_centers >= zeta_min) & \
                    (self.heat_map.zeta_centers <= zeta_max)
        region_flux = self.heat_map.flux_grid[np.ix_(theta_mask, zeta_mask)]
        return float(np.max(region_flux)) if region_flux.size > 0 else 0.0

    def _port_overlaps_exclusion(
        self, theta_c: float, zeta_c: float, theta_w: float, zeta_w: float
    ) -> bool:
        """Check if port rectangle overlaps any exclusion zone."""
        # Check corners and center of the port
        half_th = theta_w / 2
        half_ze = zeta_w / 2
        test_points = [
            (theta_c, zeta_c),                      # center
            (theta_c - half_th, zeta_c - half_ze),  # bottom-left
            (theta_c - half_th, zeta_c + half_ze),  # bottom-right
            (theta_c + half_th, zeta_c - half_ze),  # top-left
            (theta_c + half_th, zeta_c + half_ze),  # top-right
        ]
        for theta, zeta in test_points:
            if self._in_exclusion_zone(theta, zeta):
                return True
        return False

    def _in_exclusion_zone(self, theta: float, zeta: float) -> bool:
        """Check if point is in any exclusion zone."""
        for t_min, t_max, z_min, z_max in self._exclusion_zones:
            if t_min <= theta <= t_max and z_min <= zeta <= z_max:
                return True
        return False

    def solve(self, method: str = "differential_evolution") -> OptimizationResult:
        """
        Solve the port placement optimization.

        Args:
            method: Optimization method (differential_evolution or minimize)

        Returns:
            OptimizationResult with optimized port positions
        """
        try:
            from scipy.optimize import differential_evolution, minimize
        except ImportError:
            raise ImportError("scipy required for optimization: pip install scipy")

        if not self._ports:
            return OptimizationResult(
                ports=[],
                total_flux_on_ports=0.0,
                max_flux_on_ports=0.0,
                success=True,
                message="No ports to optimize",
            )

        bounds, x0 = self._build_optimization_bounds()
        x_opt, success, message = self._run_optimizer(method, bounds, x0)
        return self._build_result(x_opt, success, message)

    def _build_optimization_bounds(self) -> tuple[list, list]:
        """Build bounds and initial guess for optimizer."""
        bounds = []
        x0 = []

        for port_spec in self._ports:
            theta_w = port_spec["theta_width"]
            zeta_w = port_spec["zeta_width"]
            theta_init = port_spec["theta_init"] or 0.0
            zeta_init = port_spec["zeta_init"] or np.pi

            bounds.append((-np.pi + theta_w / 2, np.pi - theta_w / 2))
            bounds.append((zeta_w / 2, 2 * np.pi - zeta_w / 2))
            x0.extend([theta_init, zeta_init])

        return bounds, x0

    def _run_optimizer(
        self, method: str, bounds: list, x0: list
    ) -> tuple[np.ndarray, bool, str]:
        """Run the selected optimization method."""
        from scipy.optimize import differential_evolution, minimize

        if method == "differential_evolution":
            result = differential_evolution(
                self._objective,
                bounds=bounds,
                seed=42,
                maxiter=100,
                tol=1e-4,
                polish=True,
            )
        else:
            result = minimize(
                self._objective,
                x0=np.array(x0),
                bounds=bounds,
                method="L-BFGS-B",
            )

        return result.x, result.success, str(result.message)

    def _build_result(
        self, x_opt: np.ndarray, success: bool, message: str
    ) -> OptimizationResult:
        """Build OptimizationResult from optimizer output."""
        ports = []
        max_flux = 0.0
        total_flux = 0.0

        for i, port_spec in enumerate(self._ports):
            port = Port(
                name=port_spec["name"],
                theta_center=x_opt[2 * i],
                zeta_center=x_opt[2 * i + 1],
                theta_width=port_spec["theta_width"],
                zeta_width=port_spec["zeta_width"],
            )
            ports.append(port)

            flux = self.heat_map.integrated_flux_in_region(
                port.theta_min, port.theta_max,
                port.zeta_min, port.zeta_max,
            )
            total_flux += flux

            local_max = self._get_peak_flux_in_region(
                port.theta_min, port.theta_max,
                port.zeta_min, port.zeta_max,
            )
            max_flux = max(max_flux, local_max)

        return OptimizationResult(
            ports=ports,
            total_flux_on_ports=total_flux,
            max_flux_on_ports=max_flux,
            success=success,
            message=message,
        )


def plot_heat_flux_2d(
    heat_map: WallHeatMap,
    ax=None,
    cmap: str = "hot",
    vmax: Optional[float] = None,
    show_colorbar: bool = True,
    ports: Optional[list[Port]] = None,
):
    """
    Plot 2D heat flux distribution in (theta, zeta) coordinates.

    Args:
        heat_map: WallHeatMap to plot
        ax: Matplotlib axes (created if None)
        cmap: Colormap name
        vmax: Maximum value for colorbar
        show_colorbar: Whether to show colorbar
        ports: Optional list of ports to overlay

    Returns:
        Matplotlib axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))

    if vmax is None:
        vmax = heat_map.peak_flux

    im = ax.pcolormesh(
        np.degrees(heat_map.zeta_edges),
        np.degrees(heat_map.theta_edges),
        heat_map.flux_grid,
        cmap=cmap,
        vmin=0,
        vmax=vmax,
        shading="flat",
    )

    if show_colorbar:
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Heat flux (MW/m$^2$)")

    ax.set_xlabel("Toroidal angle (degrees)")
    ax.set_ylabel("Poloidal angle (degrees)")
    ax.set_title(
        f"Wall Heat Flux: {heat_map.n_lost} lost, "
        f"Peak = {heat_map.peak_flux:.2f} MW/m$^2$"
    )

    if ports:
        for port in ports:
            rect = plt.Rectangle(
                (np.degrees(port.zeta_min), np.degrees(port.theta_min)),
                np.degrees(port.zeta_width),
                np.degrees(port.theta_width),
                fill=False,
                edgecolor="cyan",
                linewidth=2,
                linestyle="--",
            )
            ax.add_patch(rect)
            ax.text(
                np.degrees(port.zeta_center),
                np.degrees(port.theta_max) + 5,
                port.name,
                ha="center",
                va="bottom",
                color="cyan",
                fontsize=10,
            )

    return ax


def plot_heat_flux_with_coil_clearance(
    heat_map: WallHeatMap,
    *,
    coil_file: str | Path,
    chartmap_file: str | Path,
    min_clearance_m: float,
    ax=None,
    cmap: str = "hot",
    vmax: Optional[float] = None,
    show_colorbar: bool = True,
    ports: Optional[list[Port]] = None,
    forbidden_alpha: float = 0.35,
):
    """
    Plot heat flux with an overlaid forbidden region based on coil clearance.

    The forbidden region is where the minimum distance from the wall to any
    coil segment is less than min_clearance_m.
    """
    import matplotlib.pyplot as plt

    if min_clearance_m <= 0.0:
        raise ValueError("min_clearance_m must be positive")

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))

    ax = plot_heat_flux_2d(
        heat_map,
        ax=ax,
        cmap=cmap,
        vmax=vmax,
        show_colorbar=show_colorbar,
        ports=ports,
    )

    wall = WallSurface.from_chartmap(chartmap_file)
    coils = CoilSegments.from_file(coil_file)
    clearance = clearance_map_on_heatmap_grid(
        wall,
        coils,
        heat_map.theta_centers,
        heat_map.zeta_centers,
    )
    forbidden = clearance < min_clearance_m

    forbidden_grid = forbidden.astype(float)
    ax.pcolormesh(
        np.degrees(heat_map.zeta_edges),
        np.degrees(heat_map.theta_edges),
        forbidden_grid,
        cmap="Greys",
        vmin=0.0,
        vmax=1.0,
        shading="flat",
        alpha=forbidden_alpha,
    )
    ax.set_title(
        f"Wall Heat Flux + Coil Clearance (d > {min_clearance_m:.3f} m)"
    )
    return ax


def plot_heat_flux_3d(
    heat_map: WallHeatMap,
    chartmap_file: Optional[str | Path] = None,
    cmap: str = "hot",
    clim: Optional[tuple[float, float]] = None,
    show: bool = True,
):
    """
    Plot 3D heat flux on wall surface using PyVista.

    Args:
        heat_map: WallHeatMap to plot
        chartmap_file: Optional chartmap NetCDF for exact wall geometry
        cmap: Colormap name
        clim: Color limits (min, max)
        show: Whether to display the plot

    Returns:
        PyVista plotter object
    """
    try:
        import pyvista as pv
    except ImportError:
        raise ImportError("pyvista required for 3D plots: pip install pyvista")

    R0 = heat_map.major_radius * 1e-2
    a = heat_map.minor_radius * 1e-2

    theta = heat_map.theta_centers
    zeta = heat_map.zeta_centers
    theta_grid, zeta_grid = np.meshgrid(theta, zeta, indexing='ij')

    R = R0 + a * np.cos(theta_grid)
    X = R * np.cos(zeta_grid)
    Y = R * np.sin(zeta_grid)
    Z = a * np.sin(theta_grid)

    grid = pv.StructuredGrid(X, Y, Z)
    grid["heat_flux"] = heat_map.flux_grid.flatten(order='F')

    plotter = pv.Plotter()
    plotter.add_mesh(
        grid,
        scalars="heat_flux",
        cmap=cmap,
        clim=clim,
        show_scalar_bar=True,
        scalar_bar_args={"title": "Heat Flux (MW/m^2)"},
    )
    plotter.add_axes()
    plotter.set_background("white")

    if show:
        plotter.show()

    return plotter


def export_to_vtk(heat_map: WallHeatMap, output_file: str | Path) -> Path:
    """
    Export heat flux to VTK file for ParaView visualization.

    Args:
        heat_map: WallHeatMap to export
        output_file: Output VTK file path

    Returns:
        Path to the created VTK file
    """
    try:
        import pyvista as pv
    except ImportError:
        raise ImportError("pyvista required for VTK export: pip install pyvista")

    R0 = heat_map.major_radius * 1e-2
    a = heat_map.minor_radius * 1e-2

    theta = heat_map.theta_centers
    zeta = heat_map.zeta_centers
    theta_grid, zeta_grid = np.meshgrid(theta, zeta, indexing='ij')

    R = R0 + a * np.cos(theta_grid)
    X = R * np.cos(zeta_grid)
    Y = R * np.sin(zeta_grid)
    Z = a * np.sin(theta_grid)

    grid = pv.StructuredGrid(X, Y, Z)
    grid["heat_flux"] = heat_map.flux_grid.flatten(order='F')
    grid["hit_count"] = heat_map.hit_count.flatten(order='F')

    output_path = Path(output_file)
    grid.save(str(output_path))
    return output_path
