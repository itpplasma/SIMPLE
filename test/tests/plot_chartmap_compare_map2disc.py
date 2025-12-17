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

    r = np.sqrt(x**2 + y**2)
    return rho, theta, zeta, r, z


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


def _closed(arr: np.ndarray) -> np.ndarray:
    return np.concatenate([arr, arr[:1]])


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="plot_chartmap_compare_map2disc")
    p.add_argument("--vmec", type=Path, required=True, help="VMEC->chartmap NetCDF file")
    p.add_argument("--map2disc", type=Path, required=True, help="map2disc chartmap NetCDF file")
    p.add_argument("--out", type=Path, required=True, help="Output PNG path")
    p.add_argument(
        "--rho",
        type=float,
        nargs="*",
        default=[0.2, 0.4, 0.6, 0.8, 1.0],
        help="Rho levels to plot (nearest grid index is used).",
    )
    args = p.parse_args(argv)

    vmec_path = args.vmec.resolve()
    map2disc_path = args.map2disc.resolve()
    out_path = args.out.resolve()

    if not vmec_path.exists() or vmec_path.stat().st_size == 0:
        raise SystemExit(f"missing vmec chartmap: {vmec_path}")
    if not map2disc_path.exists() or map2disc_path.stat().st_size == 0:
        raise SystemExit(f"missing map2disc chartmap: {map2disc_path}")

    rho_v, theta_v, zeta_v, r_v, z_v = _load_boundary_rz(vmec_path)
    rho_m, theta_m, zeta_m, r_m, z_m = _load_boundary_rz(map2disc_path)

    if rho_v.shape != rho_m.shape or not np.allclose(rho_v, rho_m, rtol=0.0, atol=0.0):
        raise SystemExit("rho grids differ between chartmaps")
    if theta_v.shape != theta_m.shape or not np.allclose(theta_v, theta_m, rtol=0.0, atol=0.0):
        raise SystemExit("theta grids differ between chartmaps")
    if zeta_v.shape != zeta_m.shape or not np.allclose(zeta_v, zeta_m, rtol=0.0, atol=0.0):
        raise SystemExit("zeta grids differ between chartmaps")

    iz_list = _select_zeta_indices(int(zeta_v.size))

    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise SystemExit(f"matplotlib required for plot test: {exc}") from exc

    rho_levels = [float(x) for x in args.rho]
    rho_levels = [x for x in rho_levels if 0.0 <= x <= 1.0]
    if not rho_levels:
        raise SystemExit("no valid --rho levels in [0,1]")

    ir_list: list[int] = []
    for rho_target in rho_levels:
        ir = int(np.argmin(np.abs(rho_v - rho_target)))
        if ir not in ir_list:
            ir_list.append(ir)

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

        for ir in ir_list:
            rv = _closed(r_v[iz, :, ir])
            zv = _closed(z_v[iz, :, ir])
            rm = _closed(r_m[iz, :, ir])
            zm = _closed(z_m[iz, :, ir])

            lab_v = "VMEC->chartmap" if ir == ir_list[-1] else None
            lab_m = "map2disc" if ir == ir_list[-1] else None

            ax_rz.plot(rv, zv, "-", lw=1.0, alpha=0.85, label=lab_v)
            ax_rz.plot(rm, zm, "--", lw=1.0, alpha=0.85, label=lab_m)
        ax_rz.set_aspect("equal", adjustable="box")
        ax_rz.set_xlabel("R [cm]")
        ax_rz.set_ylabel("Z [cm]")
        ax_rz.set_title(f"zeta={zeta_v[iz]:.6f} rad (multiple rho contours)")
        if row == 0:
            ax_rz.legend(loc="best", fontsize=9)

        ir = ir_list[-1]
        rv = r_v[iz, :, ir]
        zv = z_v[iz, :, ir]
        rm = r_m[iz, :, ir]
        zm = z_m[iz, :, ir]

        dr = rm - rv
        dz = zm - zv
        ax_d.plot(theta_v, dr, "-", lw=1.0, label="ΔR (map2disc - VMEC)")
        ax_d.plot(theta_v, dz, "-", lw=1.0, label="ΔZ (map2disc - VMEC)")
        ax_d.axhline(0.0, color="k", lw=0.6)
        ax_d.set_xlabel("theta [rad]")
        ax_d.set_ylabel("Δ [cm]")
        ax_d.set_title(
            "rho≈{:.3f}: max|ΔR|={:.3e} cm  max|ΔZ|={:.3e} cm".format(
                float(rho_v[ir]),
                float(np.max(np.abs(dr))),
                float(np.max(np.abs(dz))),
            )
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
