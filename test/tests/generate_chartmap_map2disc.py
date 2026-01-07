#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def _force_rho_convention_rho_tor(nc_path: Path) -> None:
    try:
        from netCDF4 import Dataset
    except Exception as exc:  # pragma: no cover
        raise SystemExit(f"netCDF4 required to patch rho_convention: {exc}") from exc

    with Dataset(nc_path, "a") as ds:
        current = getattr(ds, "rho_convention", "")
        if str(current).strip().lower() != "rho_tor":
            ds.setncattr("rho_convention", "rho_tor")


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
    p.add_argument(
        "--boundary-param",
        type=str,
        default="auto",
        choices=("auto", "arc", "theta"),
        help="Boundary parameterization for map2disc (auto, arc, or theta).",
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
    import inspect
    import sys

    sys.path.insert(0, str(libneo_python))

    from libneo.chartmap import write_chartmap_from_vmec_boundary

    out.parent.mkdir(parents=True, exist_ok=True)
    kwargs: dict[str, object] = {
        "nrho": int(args.nrho),
        "ntheta": int(args.ntheta),
        "nzeta": int(args.nzeta),
        "s_boundary": float(args.s_boundary),
        "M": int(args.M),
        "Nt": int(args.Nt),
        "Ng": (int(args.Ng[0]), int(args.Ng[1])),
    }
    sig_params = inspect.signature(write_chartmap_from_vmec_boundary).parameters

    boundary_offset = float(args.boundary_offset)
    if "boundary_param" in sig_params:
        boundary_param = str(args.boundary_param).strip().lower()
        if boundary_param == "auto":
            boundary_param = "theta" if boundary_offset == 0.0 else "arc"
        kwargs["boundary_param"] = boundary_param
    if boundary_offset != 0.0:
        if "boundary_offset" not in sig_params:
            raise SystemExit(
                "libneo.chartmap.write_chartmap_from_vmec_boundary does not support "
                "boundary_offset; update libneo"
            )
        kwargs["boundary_offset"] = boundary_offset
    write_chartmap_from_vmec_boundary(wout, out, **kwargs)

    _force_rho_convention_rho_tor(out)

    if not out.exists() or out.stat().st_size == 0:
        raise SystemExit(f"failed to create: {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
