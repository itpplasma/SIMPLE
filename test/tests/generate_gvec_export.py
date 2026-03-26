from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a SIMPLE GVEC export from a VMEC wout file."
    )
    parser.add_argument("--gvec-bin", required=True, type=Path)
    parser.add_argument("--exporter", required=True, type=Path)
    parser.add_argument("--wout", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--project-name", default="simple_gvec_export")
    parser.add_argument("--family", choices=["logical", "boozer"], default="logical")
    parser.add_argument("--ns", type=int, default=63)
    parser.add_argument("--ntheta", type=int, default=65)
    parser.add_argument("--nvarphi", type=int, default=65)
    parser.add_argument("--mnfactor", type=int, default=1)
    return parser.parse_args()


def write_param_file(path: Path, project_name: str, wout: Path) -> None:
    lines = [
        "! GVEC test export from VMEC",
        "whichInitEquilibrium = 1",
        f"VMECwoutfile = {wout}",
        "VMECwoutfile_format = 0",
        "sgrid_nelems = 21",
        "sgrid_grid_type = 4",
        "sgrid_rmin = 1.0e-6",
        "sgrid_rmax = 0.99",
        "X1X2_deg = 5",
        "LA_deg = 5",
        f"ProjectName = {project_name}",
        "outputIter = 0",
        "maxIter = 0",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    workdir = args.output.parent
    workdir.mkdir(parents=True, exist_ok=True)

    param_file = workdir / f"{args.project_name}.ini"
    write_param_file(param_file, args.project_name, args.wout.resolve())

    subprocess.run(
        [str(args.gvec_bin.resolve()), str(param_file.name)],
        cwd=workdir,
        check=True,
    )

    state_file = workdir / f"{args.project_name}_State_0000_00000000.dat"
    subprocess.run(
        [
            sys.executable,
            str(args.exporter.resolve()),
            "--param-file",
            str(param_file),
            "--state-file",
            str(state_file),
            "--output",
            str(args.output),
            "--family",
            args.family,
            "--ns",
            str(args.ns),
            "--ntheta",
            str(args.ntheta),
            "--nvarphi",
            str(args.nvarphi),
            "--mnfactor",
            str(args.mnfactor),
        ],
        cwd=workdir,
        check=True,
    )


if __name__ == "__main__":
    main()
