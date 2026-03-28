#!/usr/bin/env python3
"""Smoke-test the figure-8 GVEC chartmap path and generate a review plot."""

import os
import shutil
import subprocess

import netCDF4
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BUILD_DIR = os.path.join(os.path.dirname(os.path.dirname(SCRIPT_DIR)), "build")
SIMPLE_X = os.path.join(BUILD_DIR, "simple.x")
TEST_DATA = os.path.join(os.path.dirname(SCRIPT_DIR), "test_data")
CHARTMAP = os.path.join(TEST_DATA, "figure8.gvec.chartmap.nc")
BOUNDARY = os.path.join(TEST_DATA, "figure8.quasr.boundary.nc")


def write_inputs(run_dir):
    with open(os.path.join(run_dir, "start.dat"), "w") as handle:
        handle.write("0.30 0.00 0.00 1.0 0.40\n")
        handle.write("0.30 2.10 0.00 1.0 0.20\n")
        handle.write("0.45 4.20 0.10 1.0 -0.30\n")

    with open(os.path.join(run_dir, "simple.in"), "w") as handle:
        handle.write(
            """&config
multharm = 5
contr_pp = -1e10
trace_time = 1d-4
ntestpart = 3
field_input = 'figure8.gvec.chartmap.nc'
coord_input = 'figure8.gvec.chartmap.nc'
isw_field_type = 2
deterministic = .True.
integmode = 1
startmode = 2
facE_al = 1.0d0
/
"""
        )


def run_simple(run_dir):
    result = subprocess.run(
        [SIMPLE_X],
        cwd=run_dir,
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    if result.returncode != 0:
        print(result.stdout[-4000:])
        print(result.stderr[-4000:])
        raise SystemExit(result.returncode)


def load_chartmap_boundary(path):
    with netCDF4.Dataset(path) as ds:
        x = ds["x"][:].transpose(2, 1, 0)
        y = ds["y"][:].transpose(2, 1, 0)
        z = ds["z"][:].transpose(2, 1, 0)
        nfp = int(ds["num_field_periods"][:])
        return (
            x[-1] * 1e-2,
            y[-1] * 1e-2,
            z[-1] * 1e-2,
            x[0] * 1e-2,
            y[0] * 1e-2,
            z[0] * 1e-2,
            nfp,
            ds["Bmod"][:],
        )


def load_quasr_boundary(path):
    with netCDF4.Dataset(path) as ds:
        pos = ds["pos"][:]
    return pos[:, :, 0], pos[:, :, 1], pos[:, :, 2]


def tile_field_periods(xsurf, ysurf, zsurf, nfp):
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


def plot_review(run_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xb_true, yb_true, zb_true = load_quasr_boundary(BOUNDARY)
    xb, yb, zb, xa, ya, za, nfp, bmod = load_chartmap_boundary(CHARTMAP)
    xb, yb, zb = tile_field_periods(xb, yb, zb, nfp)
    xa, ya, za = tile_field_periods(xa, ya, za, nfp)
    axis_true = np.stack(
        [xb_true.mean(axis=0), yb_true.mean(axis=0), zb_true.mean(axis=0)],
        axis=1,
    )
    axis_chart = np.stack([xa.mean(axis=0), ya.mean(axis=0), za.mean(axis=0)], axis=1)

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))

    def plot_surface_projection(ax, xsurf, ysurf, zsurf, use_y, xlabel, ylabel, title):
        for j in range(0, xsurf.shape[1], max(1, xsurf.shape[1] // 10)):
            ax.plot(xsurf[:, j], ysurf[:, j] if use_y else zsurf[:, j], color="C0",
                    lw=0.8, alpha=0.6)
        for i in range(0, xsurf.shape[0], max(1, xsurf.shape[0] // 12)):
            yproj = ysurf[i] if use_y else zsurf[i]
            ax.plot(xsurf[i], yproj, color="C3", lw=0.6, alpha=0.4)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    plot_surface_projection(axes[0, 0], xb_true, yb_true, zb_true, False, "x [m]", "z [m]",
                            "QUASR boundary x-z")
    plot_surface_projection(axes[0, 1], yb_true, xb_true, zb_true, False, "y [m]", "z [m]",
                            "QUASR boundary y-z")
    plot_surface_projection(axes[0, 2], xb_true, yb_true, zb_true, True, "x [m]", "y [m]",
                            "QUASR boundary x-y")

    plot_surface_projection(axes[1, 0], xb, yb, zb, False, "x [m]", "z [m]",
                            "Chartmap surface x-z")
    axes[1, 0].plot(axis_true[:, 0], axis_true[:, 2], color="black", lw=1.6,
                    label="QUASR centerline")
    axes[1, 0].plot(axis_chart[:, 0], axis_chart[:, 2], color="goldenrod", lw=1.2,
                    ls="--", label="Chartmap rho~0")
    axes[1, 0].legend(loc="best", fontsize=8)

    plot_surface_projection(axes[1, 1], yb, xb, zb, False, "y [m]", "z [m]",
                            "Chartmap surface y-z")
    axes[1, 1].plot(axis_true[:, 1], axis_true[:, 2], color="black", lw=1.6)
    axes[1, 1].plot(axis_chart[:, 1], axis_chart[:, 2], color="goldenrod", lw=1.2, ls="--")

    ax = axes[1, 2]
    mid = bmod.shape[0] // 2
    im = ax.imshow(bmod[mid, :, :].T, origin="lower", aspect="auto", cmap="magma")
    ax.set_title("|B| at mid zeta")
    ax.set_xlabel("zeta index")
    ax.set_ylabel("theta index")
    fig.colorbar(im, ax=ax, shrink=0.8)

    fig.tight_layout()
    outpath = os.path.join(run_dir, "figure8_review.png")
    fig.savefig(outpath, dpi=180)
    print(f"Saved {outpath}")


def main():
    if not os.path.exists(SIMPLE_X):
        raise SystemExit(f"Missing simple.x: {SIMPLE_X}")
    for path in (CHARTMAP, BOUNDARY):
        if not os.path.exists(path):
            raise SystemExit(f"Missing fixture: {path}")

    run_dir = "/tmp/figure8_chartmap_smoke"
    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(run_dir)

    shutil.copy(CHARTMAP, os.path.join(run_dir, "figure8.gvec.chartmap.nc"))
    write_inputs(run_dir)
    run_simple(run_dir)

    data = np.loadtxt(os.path.join(run_dir, "confined_fraction.dat"))
    if data.ndim == 1:
        data = data[np.newaxis, :]
    confined = data[:, 1:3]
    if not np.all(np.isfinite(confined)):
        raise SystemExit("Non-finite confined fractions in figure-8 smoke test")
    if np.any(confined < -1e-12) or np.any(confined > 1.0 + 1e-12):
        raise SystemExit("Confined fractions outside [0, 1] in figure-8 smoke test")

    plot_review(run_dir)
    print("FIGURE-8 CHARTMAP PASS")


if __name__ == "__main__":
    main()
