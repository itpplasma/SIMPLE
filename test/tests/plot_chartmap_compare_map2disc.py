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


def _sort_theta(theta: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    theta = np.asarray(theta, dtype=float)
    values = np.asarray(values, dtype=float)
    theta = np.mod(theta, 2.0 * np.pi)
    order = np.argsort(theta)
    theta_s = theta[order]
    values_s = values[order]

    # Guard against duplicate theta entries (keep first occurrence).
    theta_u, idx = np.unique(theta_s, return_index=True)
    return theta_u, values_s[idx]


def _interp_periodic(theta_src: np.ndarray, values_src: np.ndarray, theta_tgt: np.ndarray) -> np.ndarray:
    theta_src = np.asarray(theta_src, dtype=float)
    values_src = np.asarray(values_src, dtype=float)
    theta_tgt = np.asarray(theta_tgt, dtype=float)

    if theta_src.size < 2:
        return np.full_like(theta_tgt, float(values_src[0]))

    theta_src_s, values_src_s = _sort_theta(theta_src, values_src)
    theta_tgt_m = np.mod(theta_tgt, 2.0 * np.pi)

    theta_ext = np.concatenate([theta_src_s, theta_src_s + 2.0 * np.pi])
    values_ext = np.concatenate([values_src_s, values_src_s])
    return np.interp(theta_tgt_m, theta_ext, values_ext)


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

    # Grids may differ slightly between generators (theta or zeta sampling).
    # For visualization we only require that the files are internally consistent.

    iz_list = _select_zeta_indices(int(zeta_v.size))

    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise SystemExit(f"matplotlib required for plot test: {exc}") from exc

    rho_levels = [float(x) for x in args.rho]
    rho_levels = [x for x in rho_levels if 0.0 <= x <= 1.0]
    if not rho_levels:
        raise SystemExit("no valid --rho levels in [0,1]")

    ir_list_v: list[int] = []
    ir_list_m: list[int] = []
    for rho_target in rho_levels:
        ir_v = int(np.argmin(np.abs(rho_v - rho_target)))
        ir_m = int(np.argmin(np.abs(rho_m - rho_target)))
        if ir_v not in ir_list_v:
            ir_list_v.append(ir_v)
        if ir_m not in ir_list_m:
            ir_list_m.append(ir_m)

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

        zeta_target = float(zeta_v[iz])
        iz_m = int(np.argmin(np.abs(zeta_m - zeta_target)))

        for ir_v, ir_m in zip(ir_list_v, ir_list_m, strict=False):
            rv = _closed(r_v[iz, :, ir_v])
            zv = _closed(z_v[iz, :, ir_v])
            rm = _closed(r_m[iz_m, :, ir_m])
            zm = _closed(z_m[iz_m, :, ir_m])

            lab_v = "VMEC->chartmap" if ir_v == ir_list_v[-1] else None
            lab_m = "map2disc" if ir_m == ir_list_m[-1] else None

            ax_rz.plot(rv, zv, "-", lw=1.0, alpha=0.85, label=lab_v)
            ax_rz.plot(rm, zm, "--", lw=1.0, alpha=0.85, label=lab_m)
        ax_rz.set_aspect("equal", adjustable="box")
        ax_rz.set_xlabel("R [cm]")
        ax_rz.set_ylabel("Z [cm]")
        ax_rz.set_title(f"zeta={zeta_target:.6f} rad (multiple rho contours)")
        if row == 0:
            ax_rz.legend(loc="best", fontsize=9)

        ir_v = ir_list_v[-1]
        ir_m = ir_list_m[-1]
        theta_target = theta_v

        rv = r_v[iz, :, ir_v]
        zv = z_v[iz, :, ir_v]
        rm = _interp_periodic(theta_m, r_m[iz_m, :, ir_m], theta_target)
        zm = _interp_periodic(theta_m, z_m[iz_m, :, ir_m], theta_target)

        dr = rm - rv
        dz = zm - zv
        ax_d.plot(theta_v, dr, "-", lw=1.0, label="ΔR (map2disc - VMEC)")
        ax_d.plot(theta_v, dz, "-", lw=1.0, label="ΔZ (map2disc - VMEC)")
        ax_d.axhline(0.0, color="k", lw=0.6)
        ax_d.set_xlabel("theta [rad]")
        ax_d.set_ylabel("Δ [cm]")
        ax_d.set_title(
            "rho≈{:.3f}: max|ΔR|={:.3e} cm  max|ΔZ|={:.3e} cm".format(
                float(rho_v[ir_v]),
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
