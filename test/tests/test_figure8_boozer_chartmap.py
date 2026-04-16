#!/usr/bin/env python3
"""Figure-8 QUASR -> GVEC -> Boozer chartmap end-to-end benchmark."""

from __future__ import annotations

import os
from pathlib import Path
import shutil

import netCDF4
import numpy as np

from boozer_chartmap_artifacts import (
    close_curve,
    close_surface_toroidally,
    confined_metrics,
    load_chartmap_fields,
    load_chartmap_surface,
    load_confined_fraction,
    plot_all_cases_total_losses,
    plot_case_loss_comparison,
    run_cmd,
    sample_indices,
    set_equal_3d_limits,
    tile_field_periods,
)


SCRIPT_DIR = Path(__file__).resolve().parent
BUILD_DIR = SCRIPT_DIR.parent.parent / "build"
SIMPLE_X = BUILD_DIR / "simple.x"
FIGURE8_DATA_ROOT = Path(
    os.environ.get(
        "SIMPLE_FIGURE8_DATA_ROOT",
        str(Path.home() / "data" / "QUASR" / "SIMPLE" / "figure8"),
    )
).expanduser()
REFERENCE_SIGNATURE = FIGURE8_DATA_ROOT / "figure8_signature_reference.txt"
BOUNDARY_FILE = FIGURE8_DATA_ROOT / "figure8.quasr.boundary.nc"
CHARTMAP_FILE = FIGURE8_DATA_ROOT / "figure8.gvec.chartmap.nc"
SIMPLE_INPUT = FIGURE8_DATA_ROOT / "simple.in"
START_INPUT = FIGURE8_DATA_ROOT / "start.dat"


def prepare_run_dir(run_dir: Path) -> tuple[Path, Path]:
    chartmap_link = run_dir / "figure8.gvec.chartmap.nc"
    chartmap_link.symlink_to(os.path.relpath(CHARTMAP_FILE, run_dir))
    shutil.copy2(SIMPLE_INPUT, run_dir / "simple.in")
    shutil.copy2(START_INPUT, run_dir / "start.dat")
    return BOUNDARY_FILE, chartmap_link


def load_quasr_boundary(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with netCDF4.Dataset(path) as ds:
        pos = np.asarray(ds["pos"][:], dtype=float)
    return pos[:, :, 0], pos[:, :, 1], pos[:, :, 2]


def figure8_signature(path: Path) -> np.ndarray:
    surface = load_chartmap_surface(path)
    fields = load_chartmap_fields(path)
    bx, by, bz = tile_field_periods(
        surface["boundary_x"], surface["boundary_y"], surface["boundary_z"], int(surface["nfp"])
    )
    ax, ay, az = tile_field_periods(
        surface["axis_x"], surface["axis_y"], surface["axis_z"], int(surface["nfp"])
    )
    axis_x = ax.mean(axis=0)
    axis_y = ay.mean(axis=0)
    axis_z = az.mean(axis=0)
    boundary_mid = bx[bx.shape[0] // 2]
    boundary_top = bz[bz.shape[0] // 4]
    idx_axis = np.linspace(0, axis_x.size - 1, 16, dtype=int)
    idx_bound = np.linspace(0, boundary_mid.size - 1, 16, dtype=int)
    rho_idx = np.linspace(0, len(fields["rho"]) - 1, 10, dtype=int)
    bmod = np.asarray(fields["Bmod"], dtype=float)
    ir = bmod.shape[0] // 2
    it = np.linspace(0, bmod.shape[1] - 1, 6, dtype=int)
    ip = np.linspace(0, bmod.shape[2] - 1, 6, dtype=int)
    sample_bmod = np.array([bmod[ir, i, j] for i in it for j in ip], dtype=float)
    return np.concatenate(
        [
            axis_x[idx_axis],
            axis_y[idx_axis],
            axis_z[idx_axis],
            boundary_mid[idx_bound],
            boundary_top[idx_bound],
            np.asarray(fields["A_phi"])[rho_idx],
            np.asarray(fields["B_theta"])[rho_idx],
            np.asarray(fields["B_phi"])[rho_idx],
            sample_bmod,
        ]
    )


def assert_figure8_golden(chartmap_path: Path, out_dir: Path) -> None:
    generated = figure8_signature(chartmap_path)
    reference = np.loadtxt(REFERENCE_SIGNATURE)
    np.savetxt(out_dir / "figure8_signature_generated.dat", generated)
    np.savetxt(out_dir / "figure8_signature_reference.dat", reference)
    if not np.allclose(generated, reference, rtol=1.0e-6, atol=1.0e-8):
        diff = np.abs(generated - reference)
        idx = int(np.argmax(diff))
        raise RuntimeError(
            "Figure-8 golden record mismatch at index "
            f"{idx}: generated={generated[idx]:.16e}, reference={reference[idx]:.16e}, "
            f"abs_diff={diff[idx]:.3e}"
        )


def plot_review(run_dir: Path, boundary_path: Path, chartmap_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xb_true, yb_true, zb_true = load_quasr_boundary(boundary_path)
    surface = load_chartmap_surface(chartmap_path)
    fields = load_chartmap_fields(chartmap_path)
    xb, yb, zb = tile_field_periods(
        surface["boundary_x"], surface["boundary_y"], surface["boundary_z"], int(surface["nfp"])
    )
    xa, ya, za = tile_field_periods(
        surface["axis_x"], surface["axis_y"], surface["axis_z"], int(surface["nfp"])
    )
    axis_true = np.stack(
        [xb_true.mean(axis=0), yb_true.mean(axis=0), zb_true.mean(axis=0)],
        axis=1,
    )
    axis_chart = np.stack([xa.mean(axis=0), ya.mean(axis=0), za.mean(axis=0)], axis=1)

    xb_true_plot, yb_true_plot, zb_true_plot = close_surface_toroidally(xb_true, yb_true, zb_true)
    xb_plot, yb_plot, zb_plot = close_surface_toroidally(xb, yb, zb)
    axis_true_plot = np.column_stack(
        [close_curve(axis_true[:, 0]), close_curve(axis_true[:, 1]), close_curve(axis_true[:, 2])]
    )
    axis_chart_plot = np.column_stack(
        [close_curve(axis_chart[:, 0]), close_curve(axis_chart[:, 1]), close_curve(axis_chart[:, 2])]
    )

    def draw_projection(ax, xsurf, ysurf, xlabel, ylabel, title, color):
        for j in sample_indices(xsurf.shape[1], max(1, xsurf.shape[1] // 12)):
            ax.plot(xsurf[:, j], ysurf[:, j], color=color, lw=0.8, alpha=0.55)
        for i in sample_indices(xsurf.shape[0], max(1, xsurf.shape[0] // 12)):
            ax.plot(xsurf[i], ysurf[i], color=color, lw=0.6, alpha=0.3)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    def draw_surface_3d(ax, xsurf, ysurf, zsurf, color, label, axis_xyz, linestyle="-"):
        for j in sample_indices(xsurf.shape[1], max(1, xsurf.shape[1] // 12)):
            ax.plot(xsurf[:, j], ysurf[:, j], zsurf[:, j], color=color, lw=0.7, alpha=0.45, ls=linestyle)
        for i in sample_indices(xsurf.shape[0], max(1, xsurf.shape[0] // 12)):
            ax.plot(xsurf[i], ysurf[i], zsurf[i], color=color, lw=0.5, alpha=0.22, ls=linestyle)
        ax.plot(axis_xyz[:, 0], axis_xyz[:, 1], axis_xyz[:, 2], color=color, lw=1.8, ls=linestyle, label=label)

    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    axes = np.empty((3, 3), dtype=object)
    for i in range(3):
        axes[0, i] = fig.add_subplot(3, 3, i + 1)
    for i in range(3):
        axes[1, i] = fig.add_subplot(3, 3, i + 4, projection="3d")
    for i in range(3):
        axes[2, i] = fig.add_subplot(3, 3, i + 7, projection="3d" if i < 2 else None)

    draw_projection(axes[0, 0], xb_true_plot, zb_true_plot, "x [m]", "z [m]", "QUASR boundary x-z", "C0")
    draw_projection(axes[0, 1], yb_true_plot, zb_true_plot, "y [m]", "z [m]", "QUASR boundary y-z", "C0")
    draw_projection(axes[0, 2], xb_true_plot, yb_true_plot, "x [m]", "y [m]", "QUASR boundary x-y", "C0")
    axes[0, 0].plot(axis_true_plot[:, 0], axis_true_plot[:, 2], color="black", lw=1.6, label="QUASR axis")
    axes[0, 0].plot(axis_chart_plot[:, 0], axis_chart_plot[:, 2], color="goldenrod", lw=1.2, ls="--", label="chartmap axis")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].plot(axis_true_plot[:, 1], axis_true_plot[:, 2], color="black", lw=1.6)
    axes[0, 1].plot(axis_chart_plot[:, 1], axis_chart_plot[:, 2], color="goldenrod", lw=1.2, ls="--")
    axes[0, 2].plot(axis_true_plot[:, 0], axis_true_plot[:, 1], color="black", lw=1.6)
    axes[0, 2].plot(axis_chart_plot[:, 0], axis_chart_plot[:, 1], color="goldenrod", lw=1.2, ls="--")

    xb_all = np.concatenate([xb_true_plot.ravel(), xb_plot.ravel()])
    yb_all = np.concatenate([yb_true_plot.ravel(), yb_plot.ravel()])
    zb_all = np.concatenate([zb_true_plot.ravel(), zb_plot.ravel()])
    views = [(22.0, -62.0), (18.0, 24.0), (78.0, -90.0)]
    for ax, (elev, azim) in zip(axes[1], views):
        draw_surface_3d(ax, xb_true_plot, yb_true_plot, zb_true_plot, "C0", "QUASR boundary", axis_true_plot)
        draw_surface_3d(ax, xb_plot, yb_plot, zb_plot, "C3", "chartmap surface", axis_chart_plot, linestyle="--")
        ax.view_init(elev=elev, azim=azim)
        set_equal_3d_limits(ax, xb_all, yb_all, zb_all)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_title(f"Overlay 3D elev={elev:.0f}, azim={azim:.0f}")
    axes[1, 0].legend(fontsize=8)

    draw_surface_3d(axes[2, 0], xb_true_plot, yb_true_plot, zb_true_plot, "C0", "QUASR boundary", axis_true_plot)
    axes[2, 0].view_init(elev=28.0, azim=-45.0)
    set_equal_3d_limits(axes[2, 0], xb_true_plot, yb_true_plot, zb_true_plot)
    axes[2, 0].set_title("QUASR boundary only")
    axes[2, 0].set_xlabel("x [m]")
    axes[2, 0].set_ylabel("y [m]")
    axes[2, 0].set_zlabel("z [m]")

    draw_surface_3d(axes[2, 1], xb_plot, yb_plot, zb_plot, "C3", "chartmap surface", axis_chart_plot)
    axes[2, 1].view_init(elev=28.0, azim=-45.0)
    set_equal_3d_limits(axes[2, 1], xb_plot, yb_plot, zb_plot)
    axes[2, 1].set_title("Chartmap surface only")
    axes[2, 1].set_xlabel("x [m]")
    axes[2, 1].set_ylabel("y [m]")
    axes[2, 1].set_zlabel("z [m]")

    bmod = np.asarray(fields["Bmod"], dtype=float)
    mid = bmod.shape[0] // 2
    im = axes[2, 2].imshow(bmod[mid].T, origin="lower", aspect="auto", cmap="magma")
    axes[2, 2].set_title("|B| at mid rho")
    axes[2, 2].set_xlabel("zeta index")
    axes[2, 2].set_ylabel("theta index")
    fig.colorbar(im, ax=axes[2, 2], shrink=0.8)

    outpath = run_dir / "figure8_review.png"
    fig.savefig(outpath, dpi=180)
    print(f"Saved {outpath}")


def maybe_plot_all_cases(run_dir: Path, figure8_metrics: dict[str, dict[str, np.ndarray | float]]) -> None:
    common_dir = Path.cwd() / "boozer_chartmap_e2e"
    if not common_dir.exists():
        return
    cases = {}
    for case_dir in sorted(common_dir.iterdir()):
        if not case_dir.is_dir():
            continue
        series = {}
        for label, subdir in [
            ("VMEC-Boozer", "vmec_run"),
            ("VMEC chartmap", "vmec_chartmap_run"),
            ("GVEC chartmap", "gvec_chartmap_run"),
        ]:
            data_file = case_dir / subdir / "confined_fraction.dat"
            if data_file.exists():
                series[label] = confined_metrics(load_confined_fraction(data_file))
        if series:
            cases[case_dir.name] = series
    cases["Figure-8"] = figure8_metrics
    plot_all_cases_total_losses(run_dir / "all_cases_losses.png", cases)


def main() -> None:
    for path, label in [
        (SIMPLE_X, "simple.x"),
        (FIGURE8_DATA_ROOT, "figure-8 data root"),
        (BOUNDARY_FILE, "figure-8 boundary"),
        (CHARTMAP_FILE, "figure-8 chartmap"),
        (SIMPLE_INPUT, "figure-8 SIMPLE input"),
        (START_INPUT, "figure-8 start file"),
        (REFERENCE_SIGNATURE, "figure-8 golden signature"),
    ]:
        if not path.exists():
            raise SystemExit(f"Missing {label}: {path}")

    run_dir = Path.cwd() / "figure8_chartmap_smoke"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)

    boundary, chartmap = prepare_run_dir(run_dir)
    run_cmd([str(SIMPLE_X)], cwd=run_dir, label="figure8 SIMPLE")

    data = load_confined_fraction(run_dir / "confined_fraction.dat")
    confined = np.asarray(data[:, 1:3], dtype=float)
    if not np.all(np.isfinite(confined)):
        raise RuntimeError("Non-finite confined fractions in figure-8 benchmark")
    if np.any(confined < -1.0e-12) or np.any(confined > 1.0 + 1.0e-12):
        raise RuntimeError("Confined fractions outside [0, 1] in figure-8 benchmark")

    metrics = {"Figure-8 GVEC chartmap": confined_metrics(data)}
    plot_case_loss_comparison(run_dir / "figure8_losses.png", "Figure-8", metrics)
    plot_review(run_dir, boundary, chartmap)
    assert_figure8_golden(chartmap, run_dir)
    maybe_plot_all_cases(run_dir, metrics)
    print("FIGURE-8 GVEC CHARTMAP PASS")


if __name__ == "__main__":
    main()
