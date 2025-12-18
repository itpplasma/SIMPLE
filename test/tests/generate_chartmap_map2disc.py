#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


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
        boundary_offset=float(args.boundary_offset),
        M=int(args.M),
        Nt=int(args.Nt),
        Ng=(int(args.Ng[0]), int(args.Ng[1])),
    )

    if not out.exists() or out.stat().st_size == 0:
        raise SystemExit(f"failed to create: {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
