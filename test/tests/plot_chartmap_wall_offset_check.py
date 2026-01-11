#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def _coord_scale_to_m(var) -> float:
    units = getattr(var, "units", "")
    units = str(units).strip().lower()
    if units in ("m", "meter", "meters"):
        return 1.0
    if units in ("cm", "centimeter", "centimeters"):
        return 0.01
    raise SystemExit(f"unsupported units '{units}' for variable '{var.name}'")


def _load_chartmap_wall_rz(
    path: Path, *, rho_target: float, zeta_indices: list[int]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with Dataset(path, "r") as ds:
        rho = np.array(ds.variables["rho"][:], dtype=float)
        theta = np.array(ds.variables["theta"][:], dtype=float)
        zeta = np.array(ds.variables["zeta"][:], dtype=float)
        x = ds.variables["x"]
        y = ds.variables["y"]
        z = ds.variables["z"]
        scale = _coord_scale_to_m(x)

        irho = int(np.argmin(np.abs(rho - rho_target)))
        r_list: list[np.ndarray] = []
        z_list: list[np.ndarray] = []
        zeta_out: list[float] = []
        for iz in zeta_indices:
            iz = int(max(0, min(zeta.size - 1, iz)))
            xx = np.array(x[iz, :, irho], dtype=float) * scale
            yy = np.array(y[iz, :, irho], dtype=float) * scale
            zz = np.array(z[iz, :, irho], dtype=float) * scale
            rr = np.sqrt(xx**2 + yy**2)
            r_list.append(rr)
            z_list.append(zz)
            zeta_out.append(float(zeta[iz]))

    r_wall = np.stack(r_list, axis=0)
    z_wall = np.stack(z_list, axis=0)
    return theta, np.array(zeta_out, dtype=float), r_wall, z_wall


def _vmec_boundary_rz(wout_path: Path, *, phi: float, ntheta: int) -> tuple[np.ndarray, np.ndarray]:
    with Dataset(wout_path, "r") as ds:
        rmnc = np.array(ds.variables["rmnc"][:], dtype=float)
        zmns = np.array(ds.variables["zmns"][:], dtype=float)
        xm = np.array(ds.variables["xm"][:], dtype=float)
        xn = np.array(ds.variables["xn"][:], dtype=float)

    js = int(rmnc.shape[0] - 1)
    theta = np.linspace(0.0, 2.0 * np.pi, int(ntheta), endpoint=False)

    angle = xm[:, None] * theta[None, :] - xn[:, None] * float(phi)
    r = (rmnc[js, :][:, None] * np.cos(angle)).sum(axis=0)
    z = (zmns[js, :][:, None] * np.sin(angle)).sum(axis=0)
    return r, z


def _select_zeta_indices(nzeta: int) -> list[int]:
    if nzeta <= 1:
        return [0]
    return list(range(int(nzeta)))


def _closed_xy(arr: np.ndarray) -> np.ndarray:
    return np.vstack([arr, arr[:1]])


def _curve_distances(
    *,
    r_wall: np.ndarray,
    z_wall: np.ndarray,
    r_vmec: np.ndarray,
    z_vmec: np.ndarray,
    offset_m: float,
) -> tuple[float, float]:
    try:
        from shapely.geometry import LinearRing, LineString, Point, Polygon
    except Exception as exc:
        raise SystemExit(f"shapely required for wall-offset test: {exc}") from exc

    wall_pts = np.column_stack([r_wall, z_wall])
    vm_pts = np.column_stack([r_vmec, z_vmec])
    poly = Polygon(LinearRing(vm_pts))
    if not poly.is_valid or poly.area <= 0.0:
        raise SystemExit("VMEC boundary polygon is invalid")

    buf = poly.buffer(float(offset_m), join_style=2, mitre_limit=2.0)
    buf_line = LineString(np.array(buf.exterior.coords, dtype=float))
    vm_line = LineString(_closed_xy(vm_pts))

    wall_to_buf = np.array([Point(p).distance(buf_line) for p in wall_pts], dtype=float)
    wall_to_vm = np.array([Point(p).distance(vm_line) for p in wall_pts], dtype=float)

    max_to_buf = float(np.max(wall_to_buf))
    max_abs_offset = float(np.max(np.abs(wall_to_vm - float(offset_m))))
    return max_to_buf, max_abs_offset


def _plot_overlay(
    *,
    out_path: Path,
    slices: list[dict[str, object]],
    offset_m: float,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"Skipping wall-offset plot overlay (matplotlib unavailable: {exc})")
        return

    try:
        n = len(slices)
        fig, axes = plt.subplots(
            nrows=1, ncols=n, figsize=(6.2 * n, 5.5), constrained_layout=True
        )
        if n == 1:
            axes = [axes]

        for ax, s in zip(axes, slices, strict=False):
            r_vm = np.asarray(s["r_vmec"], dtype=float)
            z_vm = np.asarray(s["z_vmec"], dtype=float)
            r_buf = np.asarray(s["r_buf"], dtype=float)
            z_buf = np.asarray(s["z_buf"], dtype=float)
            r_wall = np.asarray(s["r_wall"], dtype=float)
            z_wall = np.asarray(s["z_wall"], dtype=float)
            phi = float(s["phi"])

            ax.plot(r_vm, z_vm, "b-", lw=1.0, label="VMEC s=1")
            ax.plot(r_buf, z_buf, "k--", lw=1.0, label=f"buffer +{offset_m:.3f} m")
            ax.plot(r_wall, z_wall, "r.-", lw=1.0, ms=4, label="chartmap rho=1")
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("R [m]")
            ax.set_ylabel("Z [m]")
            ax.set_title(f"phi={phi:.6f} rad")
            ax.legend(loc="best", fontsize=9)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=200)
        plt.close(fig)
    except Exception as exc:
        print(f"Skipping wall-offset plot overlay (matplotlib failed: {exc})")
        return

    if not out_path.exists() or out_path.stat().st_size == 0:
        raise SystemExit(f"failed to write plot: {out_path}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="plot_chartmap_wall_offset_check")
    p.add_argument("--wout", type=Path, required=True, help="VMEC wout file")
    p.add_argument("--chartmap", type=Path, required=True, help="map2disc chartmap NetCDF file")
    p.add_argument("--out", type=Path, required=True, help="Output PNG path")
    p.add_argument("--offset", type=float, required=True, help="Expected wall offset [m]")
    p.add_argument("--rho-wall", type=float, default=1.0, help="Rho value representing wall")
    p.add_argument("--ntheta-vmec", type=int, default=4000, help="VMEC sampling resolution")
    p.add_argument("--tol-to-buffer", type=float, default=1.0e-3, help="Max dist to buffered curve [m]")
    p.add_argument("--tol-offset", type=float, default=2.0e-3, help="Max abs error in offset [m]")
    args = p.parse_args(argv)

    wout_path = args.wout.resolve()
    chartmap_path = args.chartmap.resolve()
    out_path = args.out.resolve()

    if not wout_path.exists() or wout_path.stat().st_size == 0:
        raise SystemExit(f"missing wout: {wout_path}")
    if not chartmap_path.exists() or chartmap_path.stat().st_size == 0:
        raise SystemExit(f"missing chartmap: {chartmap_path}")

    with Dataset(chartmap_path, "r") as ds:
        nzeta = int(ds.dimensions["zeta"].size)

    zeta_indices = _select_zeta_indices(nzeta)
    theta, zeta, r_wall, z_wall = _load_chartmap_wall_rz(
        chartmap_path, rho_target=float(args.rho_wall), zeta_indices=zeta_indices
    )

    worst = {
        "iz": 0,
        "phi": float(zeta[0]),
        "max_to_buffer": -1.0,
        "max_abs_offset": -1.0,
        "r_vmec": None,
        "z_vmec": None,
        "r_buf": None,
        "z_buf": None,
    }

    try:
        from shapely.geometry import LinearRing, Polygon
    except Exception as exc:
        raise SystemExit(f"shapely required for wall-offset test: {exc}") from exc

    for iz, phi in enumerate(zeta):
        r_vm, z_vm = _vmec_boundary_rz(wout_path, phi=float(phi), ntheta=int(args.ntheta_vmec))
        max_to_buffer, max_abs_offset = _curve_distances(
            r_wall=r_wall[iz],
            z_wall=z_wall[iz],
            r_vmec=r_vm,
            z_vmec=z_vm,
            offset_m=float(args.offset),
        )

        if max_to_buffer > worst["max_to_buffer"]:
            vm_pts = np.column_stack([r_vm, z_vm])
            poly = Polygon(LinearRing(vm_pts))
            buf = poly.buffer(float(args.offset), join_style=2, mitre_limit=2.0)
            bc = np.array(buf.exterior.coords, dtype=float)

            worst = {
                "iz": iz,
                "phi": float(phi),
                "max_to_buffer": float(max_to_buffer),
                "max_abs_offset": float(max_abs_offset),
                "r_vmec": r_vm,
                "z_vmec": z_vm,
                "r_buf": bc[:, 0],
                "z_buf": bc[:, 1],
            }

    if worst["max_to_buffer"] > float(args.tol_to_buffer):
        raise SystemExit(
            "chartmap wall deviates from buffered VMEC boundary: "
            f"max_to_buffer={worst['max_to_buffer']:.6e} m (tol {args.tol_to_buffer:.6e}) "
            f"at phi={worst['phi']:.6f}"
        )
    if worst["max_abs_offset"] > float(args.tol_offset):
        raise SystemExit(
            "chartmap wall offset not equal to expected value: "
            f"max_abs_offset={worst['max_abs_offset']:.6e} m (tol {args.tol_offset:.6e}) "
            f"at phi={worst['phi']:.6f}"
        )

    # Always plot phi closest to 0 and the worst-case slice (may be same).
    iz0 = int(np.argmin(np.abs(zeta - 0.0)))
    phi0 = float(zeta[iz0])
    r_vm0, z_vm0 = _vmec_boundary_rz(wout_path, phi=phi0, ntheta=int(args.ntheta_vmec))
    vm_pts0 = np.column_stack([r_vm0, z_vm0])
    poly0 = Polygon(LinearRing(vm_pts0))
    buf0 = poly0.buffer(float(args.offset), join_style=2, mitre_limit=2.0)
    bc0 = np.array(buf0.exterior.coords, dtype=float)

    slices: list[dict[str, object]] = [
        {
            "phi": phi0,
            "r_vmec": r_vm0,
            "z_vmec": z_vm0,
            "r_buf": bc0[:, 0],
            "z_buf": bc0[:, 1],
            "r_wall": r_wall[iz0],
            "z_wall": z_wall[iz0],
        }
    ]
    if int(worst["iz"]) != iz0:
        slices.append(
            {
                "phi": worst["phi"],
                "r_vmec": worst["r_vmec"],
                "z_vmec": worst["z_vmec"],
                "r_buf": worst["r_buf"],
                "z_buf": worst["z_buf"],
                "r_wall": r_wall[int(worst["iz"])],
                "z_wall": z_wall[int(worst["iz"])],
            }
        )

    _plot_overlay(out_path=out_path, slices=slices, offset_m=float(args.offset))
    print(
        "wall-offset check OK: "
        f"max_to_buffer={worst['max_to_buffer']:.3e} m, "
        f"max_abs_offset={worst['max_abs_offset']:.3e} m, "
        f"worst_phi={worst['phi']:.6f} rad"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
