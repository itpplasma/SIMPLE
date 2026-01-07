"""
Coil-aware geometric feasibility for wall ports.

This module provides:
- Parsing of SIMPLE/libneo coil point files
- Interpolation of the wall surface from a chartmap NetCDF file
- Minimum-distance queries from wall points to coil polylines

The intent is to enable port placement constraints that respect physical coil
clearances while optimizing against heat loads in (theta, zeta) coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:
    import netCDF4 as nc
    HAS_NETCDF4 = True
except ImportError:
    HAS_NETCDF4 = False


def load_coil_points(coil_file: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Load coil points from a SIMPLE/libneo-style coil file.

    Format:
        First line: integer N (number of points)
        Next N lines: x y z I

    Coordinates are expected to be in meters.

    Returns:
        points: (N, 3) array of xyz points [m]
        currents: (N,) array of currents
    """
    path = Path(coil_file)
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip()
        if not header:
            raise ValueError(f"Empty coil file: {path}")
        try:
            n_points = int(header.split()[0])
        except ValueError as exc:
            raise ValueError(f"Invalid coil file header (expected N): {path}") from exc

        data = np.loadtxt(f, dtype=float, ndmin=2)
    if data.shape[0] != n_points or data.shape[1] < 4:
        raise ValueError(
            f"Coil file {path} expected {n_points} rows of x y z I, got {data.shape}"
        )

    points = data[:, 0:3].astype(float, copy=False)
    currents = data[:, 3].astype(float, copy=False)
    return points, currents


def split_coil_polylines(
    points: np.ndarray,
    currents: np.ndarray,
    *,
    use_zero_current_separators: bool = True,
    jump_factor: float = 25.0,
) -> list[np.ndarray]:
    """
    Split concatenated coil point lists into per-coil polylines.

    Primary split rule:
        current == 0 marks a coil break (row is treated as a separator).

    Fallback split rule (for files without separators):
        Split when consecutive point spacing is a large jump compared to the
        typical segment length (jump_factor * median_length).
    """
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    if currents.ndim != 1 or currents.shape[0] != points.shape[0]:
        raise ValueError("currents must have shape (N,)")

    polylines: list[np.ndarray] = []

    if use_zero_current_separators and np.any(currents == 0.0):
        current_poly: list[np.ndarray] = []
        for p, I in zip(points, currents, strict=True):
            if I == 0.0:
                if len(current_poly) >= 2:
                    polylines.append(np.asarray(current_poly, dtype=float))
                current_poly = []
                continue
            current_poly.append(p)
        if len(current_poly) >= 2:
            polylines.append(np.asarray(current_poly, dtype=float))
        return polylines

    if points.shape[0] < 2:
        return []

    diffs = np.linalg.norm(points[1:] - points[:-1], axis=1)
    if not np.all(np.isfinite(diffs)):
        raise ValueError("Non-finite coil point coordinates encountered")

    sorted_diffs = np.sort(diffs)
    cutoff = max(1, int(0.9 * sorted_diffs.size))
    typical = float(np.median(sorted_diffs[:cutoff]))
    if typical <= 0.0:
        raise ValueError("Degenerate coil point list (zero spacing)")

    threshold = jump_factor * typical
    break_after = diffs > threshold

    start = 0
    for i, do_break in enumerate(break_after, start=0):
        if do_break:
            if i + 1 - start >= 2:
                polylines.append(points[start:i + 1].copy())
            start = i + 1
    if points.shape[0] - start >= 2:
        polylines.append(points[start:].copy())
    return polylines


@dataclass(frozen=True)
class CoilSegments:
    """A collection of polyline segments representing coils."""

    p0: np.ndarray  # (Nseg, 3)
    p1: np.ndarray  # (Nseg, 3)

    @classmethod
    def from_file(
        cls,
        coil_file: str | Path,
        *,
        use_zero_current_separators: bool = True,
        jump_factor: float = 25.0,
    ) -> "CoilSegments":
        points, currents = load_coil_points(coil_file)
        polylines = split_coil_polylines(
            points,
            currents,
            use_zero_current_separators=use_zero_current_separators,
            jump_factor=jump_factor,
        )

        if not polylines:
            raise ValueError(f"No coil polylines found in {Path(coil_file)}")

        seg_p0: list[np.ndarray] = []
        seg_p1: list[np.ndarray] = []
        for poly in polylines:
            seg_p0.append(poly[:-1])
            seg_p1.append(poly[1:])

        p0 = np.vstack(seg_p0).astype(float, copy=False)
        p1 = np.vstack(seg_p1).astype(float, copy=False)
        return cls(p0=p0, p1=p1)


def _min_distance_point_to_segments(point: np.ndarray, p0: np.ndarray, p1: np.ndarray) -> float:
    v = p1 - p0
    vv = np.einsum("ij,ij->i", v, v)
    if np.any(vv == 0.0):
        d2 = np.einsum("ij,ij->i", p0 - point, p0 - point)
        return float(np.sqrt(np.min(d2)))

    w = point - p0
    t = np.einsum("ij,ij->i", w, v) / vv
    t = np.clip(t, 0.0, 1.0)
    proj = p0 + (t[:, None] * v)
    d2 = np.einsum("ij,ij->i", proj - point, proj - point)
    return float(np.sqrt(np.min(d2)))


def min_distance_points_to_segments(
    points: np.ndarray,
    p0: np.ndarray,
    p1: np.ndarray,
    *,
    block_size: int = 128,
) -> np.ndarray:
    """
    Compute minimum distances from many points to many segments.

    Args:
        points: (M, 3) query points
        p0, p1: (Nseg, 3) segment endpoints
        block_size: segment blocking for memory control

    Returns:
        distances: (M,) minimum distance for each query point [m]
    """
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (M, 3)")
    if p0.shape != p1.shape or p0.ndim != 2 or p0.shape[1] != 3:
        raise ValueError("p0/p1 must have shape (Nseg, 3)")
    if block_size <= 0:
        raise ValueError("block_size must be positive")

    m = points.shape[0]
    nseg = p0.shape[0]
    if m == 0:
        return np.zeros((0,), dtype=float)

    if m <= 32:
        return np.asarray(
            [_min_distance_point_to_segments(points[i], p0, p1) for i in range(m)],
            dtype=float,
        )

    min_d2 = np.full((m,), np.inf, dtype=float)
    for start in range(0, nseg, block_size):
        end = min(start + block_size, nseg)
        p0b = p0[start:end]
        p1b = p1[start:end]

        v = p1b - p0b
        vv = np.einsum("ij,ij->i", v, v)
        vv = np.where(vv == 0.0, 1.0, vv)

        w = points[:, None, :] - p0b[None, :, :]
        t = np.einsum("mbi,bi->mb", w, v) / vv[None, :]
        t = np.clip(t, 0.0, 1.0)
        proj = p0b[None, :, :] + t[:, :, None] * v[None, :, :]
        d2 = np.einsum("mbi,mbi->mb", points[:, None, :] - proj, points[:, None, :] - proj)
        min_d2 = np.minimum(min_d2, np.min(d2, axis=1))

    return np.sqrt(min_d2)


def _periodic_bilinear_indices(
    grid: np.ndarray, values: np.ndarray, period: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = grid.size
    if n < 2:
        raise ValueError("Interpolation grid must have at least 2 points")

    v = np.mod(values, period)
    i1 = np.searchsorted(grid, v, side="right")
    i0 = i1 - 1
    i0 = np.mod(i0, n)
    i1 = np.mod(i1, n)

    g0 = grid[i0]
    g1 = grid[i1]

    wrapped = i1 <= i0
    denom = np.where(wrapped, (g1 + period) - g0, g1 - g0)
    denom = np.where(denom == 0.0, 1.0, denom)
    w = (v - g0) / denom
    w = np.clip(w, 0.0, 1.0)
    return i0, i1, w


@dataclass(frozen=True)
class WallSurface:
    """
    Wall surface interpolation from chartmap.

    The chartmap file stores (x, y, z) on a (rho, theta, zeta) grid for one
    field period. This class interpolates on the rho=1 surface and maps
    arbitrary zeta in [0, 2*pi) to the correct field period by rotation.
    """

    theta: np.ndarray  # (ntheta,) in [0, 2*pi)
    zeta: np.ndarray  # (nzeta,) in [0, 2*pi/nfp)
    nfp: int
    x_wall: np.ndarray  # (nzeta, ntheta) [m]
    y_wall: np.ndarray  # (nzeta, ntheta) [m]
    z_wall: np.ndarray  # (nzeta, ntheta) [m]

    @classmethod
    def from_chartmap(cls, chartmap_file: str | Path) -> "WallSurface":
        if not HAS_NETCDF4:
            raise ImportError("netCDF4 required: pip install netCDF4")

        path = Path(chartmap_file)
        with nc.Dataset(path, "r") as ds:
            theta = np.asarray(ds.variables["theta"][:], dtype=float)
            zeta = np.asarray(ds.variables["zeta"][:], dtype=float)
            nfp = int(ds.variables["num_field_periods"][...])

            x_wall = np.asarray(ds.variables["x"][:, :, -1], dtype=float) * 0.01
            y_wall = np.asarray(ds.variables["y"][:, :, -1], dtype=float) * 0.01
            z_wall = np.asarray(ds.variables["z"][:, :, -1], dtype=float) * 0.01

        return cls(
            theta=theta,
            zeta=zeta,
            nfp=nfp,
            x_wall=x_wall,
            y_wall=y_wall,
            z_wall=z_wall,
        )

    @property
    def zeta_period(self) -> float:
        return 2.0 * np.pi / float(self.nfp)

    def xyz(self, theta_hm: np.ndarray, zeta_full: np.ndarray) -> np.ndarray:
        theta_hm = np.asarray(theta_hm, dtype=float)
        zeta_full = np.asarray(zeta_full, dtype=float)

        theta_cm = np.mod(theta_hm, 2.0 * np.pi)
        period = self.zeta_period
        k = np.floor(zeta_full / period).astype(int)
        zeta_mod = zeta_full - k * period

        th0, th1, wth = _periodic_bilinear_indices(self.theta, theta_cm, 2.0 * np.pi)
        ze0, ze1, wze = _periodic_bilinear_indices(self.zeta, zeta_mod, period)

        th0, th1, wth, ze0, ze1, wze = np.broadcast_arrays(th0, th1, wth, ze0, ze1, wze)

        def interp_field(field: np.ndarray) -> np.ndarray:
            f00 = field[ze0, th0]
            f01 = field[ze0, th1]
            f10 = field[ze1, th0]
            f11 = field[ze1, th1]
            f0 = (1.0 - wth) * f00 + wth * f01
            f1 = (1.0 - wth) * f10 + wth * f11
            return (1.0 - wze) * f0 + wze * f1

        x_base = interp_field(self.x_wall)
        y_base = interp_field(self.y_wall)
        z_val = interp_field(self.z_wall)

        angle = k * period
        ca = np.cos(angle)
        sa = np.sin(angle)
        x = ca * x_base - sa * y_base
        y = sa * x_base + ca * y_base

        return np.stack([x, y, z_val], axis=-1)


@dataclass(frozen=True)
class CoilClearanceConstraint:
    wall: WallSurface
    coils: CoilSegments
    min_clearance_m: float

    def port_violates(
        self,
        theta_center: float,
        zeta_center: float,
        theta_width: float,
        zeta_width: float,
    ) -> bool:
        half_th = theta_width / 2.0
        half_ze = zeta_width / 2.0
        theta_pts = np.asarray(
            [
                theta_center,
                theta_center - half_th,
                theta_center + half_th,
                theta_center - half_th,
                theta_center + half_th,
            ],
            dtype=float,
        )
        zeta_pts = np.asarray(
            [
                zeta_center,
                zeta_center - half_ze,
                zeta_center - half_ze,
                zeta_center + half_ze,
                zeta_center + half_ze,
            ],
            dtype=float,
        )

        xyz = self.wall.xyz(theta_pts, zeta_pts).reshape((-1, 3))
        d = min_distance_points_to_segments(xyz, self.coils.p0, self.coils.p1)
        return bool(np.any(d < self.min_clearance_m))


def clearance_map_on_heatmap_grid(
    wall: WallSurface,
    coils: CoilSegments,
    theta_centers: np.ndarray,
    zeta_centers: np.ndarray,
) -> np.ndarray:
    th = np.asarray(theta_centers, dtype=float)
    ze = np.asarray(zeta_centers, dtype=float)
    th_grid, ze_grid = np.meshgrid(th, ze, indexing="ij")
    xyz = wall.xyz(th_grid, ze_grid).reshape((-1, 3))
    d = min_distance_points_to_segments(xyz, coils.p0, coils.p1)
    return d.reshape((th.size, ze.size))

