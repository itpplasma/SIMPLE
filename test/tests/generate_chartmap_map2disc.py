#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def _coord_scale_to_m(units: str) -> float:
    u = str(units).strip().lower()
    if u in ("m", "meter", "meters"):
        return 1.0
    if u in ("cm", "centimeter", "centimeters"):
        return 0.01
    raise SystemExit(f"unsupported coordinate units: {units}")


def _vmec_boundary_rz(*, wout_path: Path, phi: float, ntheta: int) -> tuple[list[float], list[float]]:
    import numpy as np
    from netCDF4 import Dataset

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
    return r.tolist(), z.tolist()


def _apply_wall_offset_inplace(*, wout: Path, chartmap: Path, offset_m: float) -> None:
    import numpy as np
    from netCDF4 import Dataset

    try:
        from shapely.geometry import LinearRing, LineString, Polygon
    except Exception as exc:
        raise SystemExit(f"shapely required for --boundary-offset: {exc}") from exc

    with Dataset(chartmap, "r+") as ds:
        rho = np.array(ds.variables["rho"][:], dtype=float)
        zeta = np.array(ds.variables["zeta"][:], dtype=float)
        x_var = ds.variables["x"]
        y_var = ds.variables["y"]
        z_var = ds.variables["z"]

        rho_wall = 1.0
        irho = int(np.argmin(np.abs(rho - float(rho_wall))))

        units = getattr(x_var, "units", "m")
        scale_units_to_m = _coord_scale_to_m(units)
        scale_m_to_units = 1.0 / scale_units_to_m

        ntheta = int(x_var.shape[1])

        for iz, phi in enumerate(zeta.tolist()):
            r_vm, z_vm = _vmec_boundary_rz(
                wout_path=wout,
                phi=float(phi),
                ntheta=max(int(ntheta), 2048),
            )
            vm_pts = np.column_stack([np.array(r_vm, dtype=float), np.array(z_vm, dtype=float)])
            poly = Polygon(LinearRing(vm_pts))
            if (not poly.is_valid) or float(poly.area) <= 0.0:
                raise SystemExit("VMEC boundary polygon is invalid; cannot apply wall offset")

            buf = poly.buffer(float(offset_m), join_style=2, mitre_limit=2.0)
            line = LineString(np.array(buf.exterior.coords, dtype=float))
            length = float(line.length)
            if length <= 0.0:
                raise SystemExit("buffered boundary has zero length; cannot apply wall offset")

            distances = np.linspace(0.0, length, num=int(ntheta) + 1, endpoint=True)[:-1]
            pts = np.array([[line.interpolate(float(d)).x, line.interpolate(float(d)).y] for d in distances])
            r_new_m = pts[:, 0]
            z_new_m = pts[:, 1]

            r_new_u = r_new_m * scale_m_to_units
            z_new_u = z_new_m * scale_m_to_units

            x_new_u = r_new_u * np.cos(float(phi))
            y_new_u = r_new_u * np.sin(float(phi))

            x_var[iz, :, irho] = x_new_u.astype(x_var.dtype, copy=False)
            y_var[iz, :, irho] = y_new_u.astype(y_var.dtype, copy=False)
            z_var[iz, :, irho] = z_new_u.astype(z_var.dtype, copy=False)


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="generate_chartmap_map2disc")
    p.add_argument(
        "--libneo-python",
        required=True,
        help="Path to libneo/python directory (for PYTHONPATH).",
    )
    p.add_argument("--wout", required=True, help="VMEC wout NetCDF path.")
    p.add_argument("--out", required=True, help="Output chartmap NetCDF path.")
    p.add_argument("--nrho", type=int, default=17)
    p.add_argument("--ntheta", type=int, default=33)
    p.add_argument("--nzeta", type=int, default=33)
    p.add_argument(
        "--s-boundary",
        type=float,
        default=1.0,
        help="VMEC surface selector s in (0, 1]; 1.0 is LCFS.",
    )
    p.add_argument(
        "--boundary-offset",
        type=float,
        default=0.0,
        help="Outward offset of VMEC boundary in meters, applied before map2disc.",
    )
    p.add_argument("--M", type=int, default=16)
    p.add_argument("--Nt", type=int, default=256)
    p.add_argument("--Ng", type=int, nargs=2, default=(256, 256))
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    libneo_python = Path(args.libneo_python).resolve()
    wout = Path(args.wout).resolve()
    out = Path(args.out).resolve()

    if not (libneo_python / "libneo").is_dir():
        raise SystemExit(f"libneo python package not found under: {libneo_python}")
    if not wout.exists():
        raise SystemExit(f"wout file not found: {wout}")

    # Import from the provided libneo source tree.
    import sys

    sys.path.insert(0, str(libneo_python))

    from libneo.chartmap import write_chartmap_from_vmec_boundary

    out.parent.mkdir(parents=True, exist_ok=True)
    write_chartmap_from_vmec_boundary(
        wout,
        out,
        nrho=int(args.nrho),
        ntheta=int(args.ntheta),
        nzeta=int(args.nzeta),
        s_boundary=float(args.s_boundary),
        M=int(args.M),
        Nt=int(args.Nt),
        Ng=(int(args.Ng[0]), int(args.Ng[1])),
    )

    if float(args.boundary_offset) != 0.0:
        _apply_wall_offset_inplace(
            wout=wout,
            chartmap=out,
            offset_m=float(args.boundary_offset),
        )

    if not out.exists() or out.stat().st_size == 0:
        raise SystemExit(f"failed to create: {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
