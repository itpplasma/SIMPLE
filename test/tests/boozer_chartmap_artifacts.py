#!/usr/bin/env python3
"""Shared artifact helpers for Boozer chartmap integration tests."""

from __future__ import annotations

from pathlib import Path
import subprocess
import urllib.request

import netCDF4
import numpy as np


STYLES = {
    "VMEC-Boozer": ("k", "-", 1.8),
    "VMEC chartmap": ("C0", "--", 1.5),
    "GVEC chartmap": ("C3", ":", 2.0),
    "Figure-8 GVEC chartmap": ("C2", "-", 1.8),
}


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def run_cmd(cmd: list[str], cwd: Path, label: str, timeout: int = 3600) -> None:
    print(f"[{label}] {' '.join(cmd)}")
    result = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    if result.returncode == 0:
        return
    print(result.stdout[-4000:])
    print(result.stderr[-4000:])
    raise RuntimeError(f"{label} failed with exit code {result.returncode}")


def download_if_missing(path: Path, url: str) -> Path:
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url} -> {path}")
    urllib.request.urlretrieve(url, path)
    return path


def load_confined_fraction(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data


def confined_metrics(data: np.ndarray) -> dict[str, np.ndarray | float]:
    time = np.asarray(data[:, 0], dtype=float)
    pass_conf = np.asarray(data[:, 1], dtype=float)
    trap_conf = np.asarray(data[:, 2], dtype=float)
    npart = float(data[0, 3])
    total_conf = pass_conf + trap_conf
    pass_lost = pass_conf[0] - pass_conf
    trap_lost = trap_conf[0] - trap_conf
    total_lost = total_conf[0] - total_conf
    return {
        "time": time,
        "npart": npart,
        "pass_conf": pass_conf,
        "trap_conf": trap_conf,
        "total_conf": total_conf,
        "pass_lost": pass_lost,
        "trap_lost": trap_lost,
        "total_lost": total_lost,
    }


def summarize_confined_fraction(data: np.ndarray) -> dict[str, float]:
    metrics = confined_metrics(data)
    return {
        "npart": float(metrics["npart"]),
        "initial_total_conf": float(metrics["total_conf"][0]),
        "final_total_conf": float(metrics["total_conf"][-1]),
        "final_total_lost": float(metrics["total_lost"][-1]),
        "final_pass_lost": float(metrics["pass_lost"][-1]),
        "final_trap_lost": float(metrics["trap_lost"][-1]),
    }


def _plot_time_axis(time: np.ndarray) -> np.ndarray:
    positive = np.asarray(time[time > 0.0], dtype=float)
    if positive.size == 0:
        return np.ones_like(time)
    floor = positive.min() / 2.0
    return np.where(time > 0.0, time, floor)


def _style(label: str) -> tuple[str, str, float]:
    return STYLES.get(label, ("0.3", "-", 1.4))


def set_equal_3d_limits(ax, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    mins = np.array([np.min(x), np.min(y), np.min(z)], dtype=float)
    maxs = np.array([np.max(x), np.max(y), np.max(z)], dtype=float)
    center = 0.5 * (mins + maxs)
    radius = 0.55 * np.max(maxs - mins)
    if radius <= 0.0:
        radius = 1.0
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.set_box_aspect((1.0, 1.0, 1.0))


def close_curve(curve: np.ndarray) -> np.ndarray:
    curve = np.asarray(curve)
    if curve.ndim != 1 or curve.size == 0:
        return curve
    return np.concatenate([curve, curve[:1]])


def close_surface_toroidally(
    xsurf: np.ndarray, ysurf: np.ndarray, zsurf: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xsurf = np.asarray(xsurf)
    ysurf = np.asarray(ysurf)
    zsurf = np.asarray(zsurf)
    if xsurf.ndim != 2 or xsurf.shape[1] == 0:
        return xsurf, ysurf, zsurf
    return (
        np.concatenate([xsurf, xsurf[:, :1]], axis=1),
        np.concatenate([ysurf, ysurf[:, :1]], axis=1),
        np.concatenate([zsurf, zsurf[:, :1]], axis=1),
    )


def sample_indices(size: int, stride: int) -> list[int]:
    indices = list(range(0, size, stride))
    last = size - 1
    if last not in indices:
        indices.append(last)
    return indices


def plot_surface_lines_3d(
    ax,
    xsurf: np.ndarray,
    ysurf: np.ndarray,
    zsurf: np.ndarray,
    color: str,
    label: str,
    axis_x: np.ndarray,
    axis_y: np.ndarray,
    axis_z: np.ndarray,
    linestyle: str = "-",
) -> None:
    xsurf, ysurf, zsurf = close_surface_toroidally(xsurf, ysurf, zsurf)
    axis_x = close_curve(axis_x)
    axis_y = close_curve(axis_y)
    axis_z = close_curve(axis_z)
    n_phi = xsurf.shape[1]
    n_rho = xsurf.shape[0]
    phi_stride = max(1, n_phi // 12)
    rho_stride = max(1, n_rho // 12)
    for idx in sample_indices(n_phi, phi_stride):
        ax.plot(
            xsurf[:, idx],
            ysurf[:, idx],
            zsurf[:, idx],
            color=color,
            lw=0.7,
            alpha=0.45,
            ls=linestyle,
        )
    for idx in sample_indices(n_rho, rho_stride):
        ax.plot(
            xsurf[idx],
            ysurf[idx],
            zsurf[idx],
            color=color,
            lw=0.5,
            alpha=0.22,
            ls=linestyle,
        )
    ax.plot(axis_x, axis_y, axis_z, color=color, lw=1.8, ls=linestyle, label=label)


def plot_case_loss_comparison(
    output: Path,
    case_name: str,
    metrics_by_label: dict[str, dict[str, np.ndarray | float]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(output)
    ref = metrics_by_label.get("VMEC-Boozer")
    if ref is None:
        ref = next(iter(metrics_by_label.values()))
    time_ms = _plot_time_axis(np.asarray(ref["time"])) * 1.0e3

    panels = [
        ("total_conf", "Total confined", False),
        ("total_lost", "Total lost", False),
        ("pass_conf", "Passing confined", False),
        ("pass_lost", "Passing lost", False),
        ("trap_conf", "Trapped confined", False),
        ("trap_lost", "Trapped lost", False),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for ax, (key, title, semilog_y) in zip(axes.flat, panels):
        for label, metrics in metrics_by_label.items():
            color, linestyle, width = _style(label)
            y = np.asarray(metrics[key], dtype=float)
            if semilog_y:
                ax.semilogy(time_ms, y, color=color, ls=linestyle, lw=width, label=label)
            else:
                ax.plot(time_ms, y, color=color, ls=linestyle, lw=width, label=label)
            ax.set_xscale("log")
        ax.set_xlabel("time [ms]")
        ax.set_ylabel("fraction")
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(title)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"{case_name}: confined and lost fractions", fontsize=14)
    fig.savefig(output, dpi=180)
    print(f"Saved {output}")


def plot_all_cases_total_losses(
    output: Path,
    cases: dict[str, dict[str, dict[str, np.ndarray | float]]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(output)
    names = list(cases)
    panels = [
        ("total_lost", "Total lost"),
        ("pass_lost", "Passing lost"),
        ("trap_lost", "Trapped lost"),
    ]
    fig, axes = plt.subplots(
        len(names),
        len(panels),
        figsize=(15, 3.4 * len(names)),
        constrained_layout=True,
    )
    if len(names) == 1:
        axes = np.asarray([axes])
    for row, case_name in enumerate(names):
        metrics_by_label = cases[case_name]
        ref = next(iter(metrics_by_label.values()))
        time_ms = _plot_time_axis(np.asarray(ref["time"])) * 1.0e3
        for col, (key, title) in enumerate(panels):
            ax = axes[row, col]
            for label, metrics in metrics_by_label.items():
                color, linestyle, width = _style(label)
                ax.plot(
                    time_ms,
                    np.asarray(metrics[key], dtype=float),
                    color=color,
                    ls=linestyle,
                    lw=width,
                    label=label,
                )
            ax.set_xscale("log")
            ax.set_ylim(-0.05, 1.05)
            ax.set_xlabel("time [ms]")
            if col == 0:
                ax.set_ylabel(case_name)
            ax.set_title(title)
            if row == 0 and col == 0:
                ax.legend(fontsize=8)
    fig.suptitle("Loss comparison across all benchmark cases", fontsize=14)
    fig.savefig(output, dpi=180)
    print(f"Saved {output}")


def load_chartmap_surface(path: Path) -> dict[str, np.ndarray | int]:
    with netCDF4.Dataset(path) as ds:
        x = np.asarray(ds["x"][:]).transpose(2, 1, 0) * 1.0e-2
        y = np.asarray(ds["y"][:]).transpose(2, 1, 0) * 1.0e-2
        z = np.asarray(ds["z"][:]).transpose(2, 1, 0) * 1.0e-2
        nfp = int(np.asarray(ds["num_field_periods"][:]).item())
    return {
        "boundary_x": x[-1],
        "boundary_y": y[-1],
        "boundary_z": z[-1],
        "axis_x": x[0],
        "axis_y": y[0],
        "axis_z": z[0],
        "nfp": nfp,
    }


def load_chartmap_fields(path: Path) -> dict[str, np.ndarray]:
    with netCDF4.Dataset(path) as ds:
        rho = np.asarray(ds["rho"][:], dtype=float)
        theta = np.asarray(ds["theta"][:], dtype=float)
        zeta = np.asarray(ds["zeta"][:], dtype=float)
        bmod = np.asarray(ds["Bmod"][:], dtype=float).transpose(2, 1, 0)
        a_phi = np.asarray(ds["A_phi"][:], dtype=float)
        b_theta = np.asarray(ds["B_theta"][:], dtype=float)
        b_phi = np.asarray(ds["B_phi"][:], dtype=float)
    return {
        "rho": rho,
        "theta": theta,
        "zeta": zeta,
        "Bmod": bmod,
        "A_phi": a_phi,
        "B_theta": b_theta,
        "B_phi": b_phi,
    }


def tile_field_periods(xsurf: np.ndarray, ysurf: np.ndarray, zsurf: np.ndarray, nfp: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_parts = []
    y_parts = []
    z_parts = []
    for k in range(nfp):
        angle = 2.0 * np.pi * k / nfp
        cang = np.cos(angle)
        sang = np.sin(angle)
        x_parts.append(cang * xsurf - sang * ysurf)
        y_parts.append(sang * xsurf + cang * ysurf)
        z_parts.append(zsurf.copy())
    return np.concatenate(x_parts, axis=1), np.concatenate(y_parts, axis=1), np.concatenate(z_parts, axis=1)


def plot_surface_comparison(
    output: Path,
    case_name: str,
    left_path: Path,
    right_path: Path,
    left_label: str,
    right_label: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(output)
    left = load_chartmap_surface(left_path)
    right = load_chartmap_surface(right_path)
    left_x, left_y, left_z = tile_field_periods(left["boundary_x"], left["boundary_y"], left["boundary_z"], int(left["nfp"]))
    right_x, right_y, right_z = tile_field_periods(right["boundary_x"], right["boundary_y"], right["boundary_z"], int(right["nfp"]))
    left_ax_x, left_ax_y, left_ax_z = tile_field_periods(left["axis_x"], left["axis_y"], left["axis_z"], int(left["nfp"]))
    right_ax_x, right_ax_y, right_ax_z = tile_field_periods(right["axis_x"], right["axis_y"], right["axis_z"], int(right["nfp"]))

    fig = plt.figure(figsize=(15, 8.5), constrained_layout=True)
    axes = np.empty((2, 3), dtype=object)
    for i in range(3):
        axes[0, i] = fig.add_subplot(2, 3, i + 1)
    for i in range(3):
        axes[1, i] = fig.add_subplot(2, 3, i + 4, projection="3d")

    panels = [
        ("x-z", left_x, left_z, right_x, right_z, left_ax_x, left_ax_z, right_ax_x, right_ax_z, "x [m]", "z [m]"),
        ("y-z", left_y, left_z, right_y, right_z, left_ax_y, left_ax_z, right_ax_y, right_ax_z, "y [m]", "z [m]"),
        ("x-y", left_x, left_y, right_x, right_y, left_ax_x, left_ax_y, right_ax_x, right_ax_y, "x [m]", "y [m]"),
    ]
    for ax, (title, xl, yl, xr, yr, xal, yal, xar, yar, xlabel, ylabel) in zip(axes[0], panels):
        xl, yl, _ = close_surface_toroidally(xl, yl, yl)
        xr, yr, _ = close_surface_toroidally(xr, yr, yr)
        for idx in sample_indices(xl.shape[1], max(1, xl.shape[1] // 12)):
            ax.plot(xl[:, idx], yl[:, idx], color="C0", lw=0.8, alpha=0.45)
        for idx in sample_indices(xr.shape[1], max(1, xr.shape[1] // 12)):
            ax.plot(xr[:, idx], yr[:, idx], color="C3", lw=0.8, alpha=0.45)
        for idx in sample_indices(xl.shape[0], max(1, xl.shape[0] // 12)):
            ax.plot(xl[idx], yl[idx], color="C0", lw=0.6, alpha=0.3)
        for idx in sample_indices(xr.shape[0], max(1, xr.shape[0] // 12)):
            ax.plot(xr[idx], yr[idx], color="C3", lw=0.6, alpha=0.3)
        ax.plot(close_curve(xal.mean(axis=0)), close_curve(yal.mean(axis=0)),
                color="C0", lw=2.0, label=left_label)
        ax.plot(close_curve(xar.mean(axis=0)), close_curve(yar.mean(axis=0)),
                color="C3", lw=1.6, ls="--", label=right_label)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    axes[0, 0].legend(fontsize=8)

    views = [(22.0, -62.0), (18.0, 22.0), (78.0, -90.0)]
    all_x = np.concatenate([left_x.ravel(), right_x.ravel(), left_ax_x.ravel(), right_ax_x.ravel()])
    all_y = np.concatenate([left_y.ravel(), right_y.ravel(), left_ax_y.ravel(), right_ax_y.ravel()])
    all_z = np.concatenate([left_z.ravel(), right_z.ravel(), left_ax_z.ravel(), right_ax_z.ravel()])
    for ax, (elev, azim) in zip(axes[1], views):
        plot_surface_lines_3d(
            ax,
            left_x,
            left_y,
            left_z,
            "C0",
            left_label,
            left_ax_x.mean(axis=0),
            left_ax_y.mean(axis=0),
            left_ax_z.mean(axis=0),
        )
        plot_surface_lines_3d(
            ax,
            right_x,
            right_y,
            right_z,
            "C3",
            right_label,
            right_ax_x.mean(axis=0),
            right_ax_y.mean(axis=0),
            right_ax_z.mean(axis=0),
            linestyle="--",
        )
        ax.view_init(elev=elev, azim=azim)
        set_equal_3d_limits(ax, all_x, all_y, all_z)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_title(f"3D view elev={elev:.0f}, azim={azim:.0f}")
    axes[1, 0].legend(fontsize=8)
    fig.suptitle(f"{case_name}: flux-surface shape comparison", fontsize=14)
    fig.savefig(output, dpi=180)
    print(f"Saved {output}")


def plot_field_component_comparison(
    output: Path,
    case_name: str,
    left_path: Path,
    right_path: Path,
    left_label: str,
    right_label: str,
) -> dict[str, float]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(output)
    left = load_chartmap_fields(left_path)
    right = load_chartmap_fields(right_path)

    rho_common = np.linspace(0.0, 1.0, min(len(left["rho"]), len(right["rho"])))
    left_a = np.interp(rho_common, left["rho"], left["A_phi"])
    right_a = np.interp(rho_common, right["rho"], right["A_phi"])
    left_bt = np.interp(rho_common, left["rho"], left["B_theta"])
    right_bt = np.interp(rho_common, right["rho"], right["B_theta"])
    left_bp = np.interp(rho_common, left["rho"], left["B_phi"])
    right_bp = np.interp(rho_common, right["rho"], right["B_phi"])

    bmod_left = np.asarray(left["Bmod"])
    bmod_right = np.asarray(right["Bmod"])
    nrho = min(bmod_left.shape[0], bmod_right.shape[0])
    ntheta = min(bmod_left.shape[1], bmod_right.shape[1])
    nzeta = min(bmod_left.shape[2], bmod_right.shape[2])
    bmod_left = bmod_left[:nrho, :ntheta, :nzeta]
    bmod_right = bmod_right[:nrho, :ntheta, :nzeta]
    mid = nrho // 2
    rel_bmod = np.abs(bmod_right[mid] - bmod_left[mid]) / np.maximum(np.abs(bmod_left[mid]), 1.0e-14)

    metrics = {
        "max_rel_A_phi": float(np.max(np.abs(right_a - left_a) / np.maximum(np.abs(left_a), 1.0e-14))),
        "max_rel_B_theta": float(np.max(np.abs(right_bt - left_bt) / np.maximum(np.abs(left_bt), 1.0e-14))),
        "max_rel_B_phi": float(np.max(np.abs(right_bp - left_bp) / np.maximum(np.abs(left_bp), 1.0e-14))),
        "max_rel_Bmod": float(np.max(rel_bmod)),
    }

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    profiles = [
        ("A_phi", left_a, right_a, "A_phi"),
        ("B_theta", left_bt, right_bt, "B_theta"),
        ("B_phi", left_bp, right_bp, "B_phi"),
    ]
    for ax, (title, left_vals, right_vals, key) in zip(axes[0], profiles):
        ax.plot(rho_common, left_vals, color="C0", lw=1.8, label=left_label)
        ax.plot(rho_common, right_vals, color="C3", lw=1.4, ls="--", label=right_label)
        ax.set_xlabel("rho")
        ax.set_title(f"{title} (max rel {metrics[f'max_rel_{key}']:.2e})")
        ax.legend(fontsize=8)

    im0 = axes[1, 0].imshow(bmod_left[mid].T, origin="lower", aspect="auto", cmap="magma")
    axes[1, 0].set_title(f"{left_label} |B|")
    axes[1, 0].set_xlabel("zeta index")
    axes[1, 0].set_ylabel("theta index")
    fig.colorbar(im0, ax=axes[1, 0], shrink=0.8)

    im1 = axes[1, 1].imshow(bmod_right[mid].T, origin="lower", aspect="auto", cmap="magma")
    axes[1, 1].set_title(f"{right_label} |B|")
    axes[1, 1].set_xlabel("zeta index")
    axes[1, 1].set_ylabel("theta index")
    fig.colorbar(im1, ax=axes[1, 1], shrink=0.8)

    im2 = axes[1, 2].imshow(np.log10(np.maximum(rel_bmod.T, 1.0e-16)), origin="lower", aspect="auto", cmap="viridis")
    axes[1, 2].set_title(f"log10 rel |B| diff (max {metrics['max_rel_Bmod']:.2e})")
    axes[1, 2].set_xlabel("zeta index")
    axes[1, 2].set_ylabel("theta index")
    fig.colorbar(im2, ax=axes[1, 2], shrink=0.8)

    fig.suptitle(f"{case_name}: chartmap field comparison", fontsize=14)
    fig.savefig(output, dpi=180)
    print(f"Saved {output}")
    return metrics
