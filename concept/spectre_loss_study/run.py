#!/usr/bin/env python3
"""Matched physical pilots through SIMPLE's native make/CTest authority."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


STUDY_SOURCES = (
    "concept/spectre_loss_study/study_spectre_orbits.f90",
    "src/simple.f90", "src/spectre_sympl_orbit.f90",
    "src/spectre_fo_hybrid.f90", "src/spectre_fo_field.f90",
    "src/orbit_symplectic.f90", "src/orbit_symplectic_base.f90",
    "src/field_can.f90", "src/field/field_can_spectre.f90",
)


def source_hashes(root):
    return {path: hashlib.sha256((root / path).read_bytes()).hexdigest()
            for path in STUDY_SOURCES if (root / path).exists()}


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--starts", type=Path)
    parser.add_argument("--field", choices=["ncsx", "g3v3"], required=True)
    parser.add_argument("--sampling", choices=["volume", "sheet", "phase"],
                        default="volume")
    parser.add_argument("--mode", choices=["hybrid", "fo"], default="hybrid")
    parser.add_argument("--grid", type=int, nargs=3)
    parser.add_argument("--markers", type=int, default=16)
    parser.add_argument("--duration", type=float, default=2e-6)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--entry", type=float, default=3)
    parser.add_argument("--samples", type=int, default=65)
    parser.add_argument("--phase", type=float, default=0,
                        help="FO launch phase offset in radians (research only)")
    parser.add_argument("--npoiper2", type=int, default=1024)
    parser.add_argument("--macrosteps", type=int, default=64)
    parser.add_argument("--construction-tolerance", type=float)
    args = parser.parse_args()
    args.root = args.root.resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    # Both revisions use the same equilibrium bytes, independently of checkout.
    source_root = Path(__file__).resolve().parents[2]
    fixture = (source_root / "test/tests/data/ncsx_freeboundary.h5"
               if args.field == "ncsx" else source_root /
               "build/_deps/libneo-src/test/spectre/data/spectre_test.h5")
    config = args.output / "study.in"
    config.write_text(
        "&config\n"
        f"field_input='{fixture}'\n"
        "integ_coords=6, integmode=3, crossing_level=3\n"
        f"ntestpart={args.markers}, ntimstep={args.macrosteps + 1}, "
        f"npoiper2={args.npoiper2}\n"
        f"facE_al=1000, trace_time={args.duration:.17e}, relerr=1d-12\n"
        f"spectre_fo_entry_larmor={1e6 if args.mode == 'fo' else args.entry}\n"
        f"spectre_fo_tolerance={args.tolerance:.17e}\n"
        + ("spectre_fo_grid=" + ",".join(map(str, args.grid)) + "\n"
           if args.grid else "")
        + (f"spectre_ncon_ode_relerr={args.construction_tolerance:.17e}\n"
           if args.construction_tolerance is not None else "")
        + "deterministic=.true., ran_seed=12345 /\n")
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1", OMP_PROC_BIND="close", OMP_PLACES="cores",
               SIMPLE_SPECTRE_STUDY_INPUT=str(config),
               SIMPLE_SPECTRE_STUDY_OUTPUT=str(args.output / "study"),
               SIMPLE_SPECTRE_STUDY_MODE=args.mode,
               SIMPLE_SPECTRE_STUDY_SAMPLES=str(args.samples),
               SIMPLE_SPECTRE_STUDY_PHASE=str(args.phase),
               SIMPLE_SPECTRE_STUDY_SAMPLING=args.sampling)
    env.pop("SIMPLE_SPECTRE_STUDY_STARTS", None)
    if args.starts:
        env["SIMPLE_SPECTRE_STUDY_STARTS"] = str(args.starts.resolve())
    command = ["make", "test", "CONFIG=Fast", "TEST=^study_spectre_orbits$",
               "VERBOSE=1"]
    started = time.monotonic()
    with (args.output / "native.log").open("w") as log:
        process = subprocess.run(command, cwd=args.root, env=env, stdout=log,
                                 stderr=subprocess.STDOUT, check=False)
    metadata = vars(args).copy()
    metadata.update(source_revision=git(args.root, "rev-parse", "HEAD"),
                    measurement_diff=git(args.root, "diff", "--stat"),
                    source_hashes=source_hashes(args.root),
                    fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),
                    process_s=time.monotonic() - started, returncode=process.returncode,
                    load_average=os.getloadavg(), command=command)
    neo_source = args.root / "build/_deps/libneo-src"
    if not (neo_source / ".git").exists():
        neo_source = source_root / "build/_deps/libneo-src"
    metadata["libneo_revision"] = git(neo_source, "rev-parse", "HEAD")
    (args.output / "run.json").write_text(json.dumps(metadata, default=str, indent=2) + "\n")
    summary = args.output / "study_summary.json"
    print(json.dumps(dict(case=str(args.output), returncode=process.returncode,
                          **(json.loads(summary.read_text()) if summary.exists() else {}))))
    raise SystemExit(process.returncode)


if __name__ == "__main__":
    main()
