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
        return (
            ds["x"][:].transpose(2, 1, 0)[-1],
            ds["y"][:].transpose(2, 1, 0)[-1],
            ds["z"][:].transpose(2, 1, 0)[-1],
            ds["Bmod"][:],
        )


def load_quasr_boundary(path):
    with netCDF4.Dataset(path) as ds:
        pos = ds["pos"][:]
    return pos[:, :, 0], pos[:, :, 1], pos[:, :, 2]


def plot_review(run_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xb_true, yb_true, zb_true = load_quasr_boundary(BOUNDARY)
    xb, yb, zb, bmod = load_chartmap_boundary(CHARTMAP)

    fig = plt.figure(figsize=(13, 4))

    ax = fig.add_subplot(1, 3, 1, projection="3d")
    for j in range(0, xb_true.shape[0], 8):
        ax.plot(xb_true[j], yb_true[j], zb_true[j], color="C0", lw=0.8, alpha=0.7)
    for i in range(0, xb_true.shape[1], 8):
        ax.plot(xb_true[:, i], yb_true[:, i], zb_true[:, i], color="C3", lw=0.6, alpha=0.45)
    ax.set_title("QUASR figure-8 boundary")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.view_init(elev=18, azim=40)

    ax = fig.add_subplot(1, 3, 2, projection="3d")
    for j in range(0, xb.shape[1], 2):
        ax.plot(xb[:, j], yb[:, j], zb[:, j], color="C0", lw=0.8, alpha=0.7)
    for i in range(0, xb.shape[0], 2):
        ax.plot(xb[i], yb[i], zb[i], color="C3", lw=0.6, alpha=0.45)
    ax.set_title("SIMPLE chartmap geometry")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.set_zlabel("z [cm]")
    ax.view_init(elev=18, azim=40)

    ax = fig.add_subplot(1, 3, 3)
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
