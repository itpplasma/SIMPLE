#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def _load_boundary_rz(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with Dataset(path, "r") as ds:
        rho = np.array(ds.variables["rho"][:], dtype=float)
        theta = np.array(ds.variables["theta"][:], dtype=float)
        zeta = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)

    ir = int(rho.size - 1)
    r = np.sqrt(x[:, :, ir] ** 2 + y[:, :, ir] ** 2)
    zz = z[:, :, ir]
    return theta, zeta, r, zz


def _select_zeta_indices(nzeta: int) -> list[int]:
    if nzeta <= 1:
        return [0]
    candidates = [0, nzeta // 4, nzeta // 2, (3 * nzeta) // 4, nzeta - 1]
    out: list[int] = []
    for i in candidates:
        ii = int(max(0, min(nzeta - 1, i)))
        if ii not in out:
            out.append(ii)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="plot_chartmap_compare_map2disc")
    p.add_argument("--vmec", type=Path, required=True, help="VMEC->chartmap NetCDF file")
    p.add_argument("--map2disc", type=Path, required=True, help="map2disc chartmap NetCDF file")
    p.add_argument("--out", type=Path, required=True, help="Output PNG path")
    args = p.parse_args(argv)

    vmec_path = args.vmec.resolve()
    map2disc_path = args.map2disc.resolve()
    out_path = args.out.resolve()

    if not vmec_path.exists() or vmec_path.stat().st_size == 0:
        raise SystemExit(f"missing vmec chartmap: {vmec_path}")
    if not map2disc_path.exists() or map2disc_path.stat().st_size == 0:
        raise SystemExit(f"missing map2disc chartmap: {map2disc_path}")

    theta_v, zeta_v, r_v, z_v = _load_boundary_rz(vmec_path)
    theta_m, zeta_m, r_m, z_m = _load_boundary_rz(map2disc_path)

    if theta_v.shape != theta_m.shape:
        raise SystemExit("theta grids differ between chartmaps")
    if zeta_v.shape != zeta_m.shape:
        raise SystemExit("zeta grids differ between chartmaps")

    iz_list = _select_zeta_indices(int(zeta_v.size))

    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise SystemExit(f"matplotlib required for plot test: {exc}") from exc

    nrows = len(iz_list)
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=2,
        figsize=(10.0, 2.4 * nrows),
        constrained_layout=True,
    )
    if nrows == 1:
        axes = np.array([axes])

    for row, iz in enumerate(iz_list):
        ax_rz = axes[row, 0]
        ax_d = axes[row, 1]

        rv = r_v[iz, :]
        zv = z_v[iz, :]
        rm = r_m[iz, :]
        zm = z_m[iz, :]

        ax_rz.plot(rv, zv, "-", lw=1.2, label="VMEC->chartmap")
        ax_rz.plot(rm, zm, "--", lw=1.2, label="map2disc")
        ax_rz.set_aspect("equal", adjustable="box")
        ax_rz.set_xlabel("R [cm]")
        ax_rz.set_ylabel("Z [cm]")
        ax_rz.set_title(f"rho=1 slice at zeta={zeta_v[iz]:.6f} rad")
        if row == 0:
            ax_rz.legend(loc="best", fontsize=9)

        dr = rm - rv
        dz = zm - zv
        ax_d.plot(theta_v, dr, "-", lw=1.0, label="ΔR (map2disc - VMEC)")
        ax_d.plot(theta_v, dz, "-", lw=1.0, label="ΔZ (map2disc - VMEC)")
        ax_d.axhline(0.0, color="k", lw=0.6)
        ax_d.set_xlabel("theta [rad]")
        ax_d.set_ylabel("Δ [cm]")
        ax_d.set_title(
            f"max|ΔR|={np.max(np.abs(dr)):.3e} cm  max|ΔZ|={np.max(np.abs(dz)):.3e} cm"
        )
        if row == 0:
            ax_d.legend(loc="best", fontsize=9)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    if not out_path.exists() or out_path.stat().st_size == 0:
        raise SystemExit(f"failed to write plot: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

