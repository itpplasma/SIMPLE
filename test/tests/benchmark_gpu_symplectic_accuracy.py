#!/usr/bin/env python3
"""Measure short fixed-particle GPU accuracy against a tight DOPRI control.

This is intentionally a small, dependency-free harness.  It compares the
particle profile emitted by SIMPLE rather than only comparing aggregate loss
counts, so a candidate cannot hide an endpoint or loss-time discrepancy.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess


def read_value(pattern: str, output: str, name: str) -> float:
    match = re.search(pattern, output)
    if match is None:
        raise RuntimeError(f"{name} missing from SIMPLE output")
    return float(match.group(1))


def read_integer(pattern: str, output: str, name: str) -> int:
    match = re.search(pattern, output)
    if match is None:
        raise RuntimeError(f"{name} missing from SIMPLE output")
    return int(match.group(1))


def angle_difference(left: float, right: float) -> float:
    return abs(math.atan2(math.sin(left - right), math.cos(left - right)))


def write_input(
    path: Path,
    trace_time: float,
    ntimstep: int,
    particle_count: int,
    npoiper2: int,
    relerr: float,
    netcdffile: Path | None,
    field_input: Path | None,
    coord_input: Path | None,
) -> None:
    lines = [
        "&config",
        f"trace_time = {trace_time:.17e}",
        f"ntimstep = {ntimstep}",
        f"ntestpart = {particle_count}",
        "startmode = 2",
        "isw_field_type = 2",
        "integmode = 1",
        f"npoiper2 = {npoiper2}",
        f"relerr = {relerr:.17e}",
    ]
    if netcdffile is not None:
        lines.append(f"netcdffile = '{netcdffile}'")
    if field_input is not None:
        lines.append(f"field_input = '{field_input}'")
    if coord_input is not None:
        lines.append(f"coord_input = '{coord_input}'")
    lines.append("/")
    path.write_text("\n".join(lines) + "\n")


def read_profile(path: Path) -> list[dict[str, float | int]]:
    with path.open(newline="") as stream:
        rows: list[dict[str, float | int]] = []
        for row in csv.DictReader(stream):
            rows.append(
                {
                    "s_end": float(row["s_end"]),
                    "theta_end": float(row["theta_end"]),
                    "zeta_end": float(row["zeta_end"]),
                    "speed_end": float(row["speed_end"]),
                    "pitch_end": float(row["pitch_end"]),
                    "loss_step": int(row["loss_step"]),
                    "loss_time": float(row["loss_time"]),
                    "nfev": int(row["nfev"]),
                }
            )
    return rows


def run_case(
    executable: Path,
    case_dir: Path,
    starts: str,
    particle_count: int,
    method: str,
    npoiper2: int,
    relerr: float,
    trace_time: float,
    ntimstep: int,
    netcdffile: Path | None,
    field_input: Path | None,
    coord_input: Path | None,
) -> dict[str, object]:
    case_dir.mkdir(parents=True)
    (case_dir / "start.dat").write_text(starts)
    write_input(
        case_dir / "simple.in",
        trace_time,
        ntimstep,
        particle_count,
        npoiper2,
        relerr,
        netcdffile,
        field_input,
        coord_input,
    )
    profile_path = case_dir / "particles.csv"
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "SIMPLE_GPU_RUN": "1",
            "SIMPLE_GPU_METHOD": method,
            "SIMPLE_GPU_START_COORDINATES": "boozer",
            "SIMPLE_GPU_LANDREMAN": "0",
            "SIMPLE_GPU_PARTICLE_PROFILE": str(profile_path),
        }
    )
    result = subprocess.run(
        [executable, "simple.in"],
        cwd=case_dir,
        env=environment,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (case_dir / "run.log").write_text(result.stdout)
    return {
        "method": method,
        "npoiper2": npoiper2,
        "relerr": relerr,
        "tracing_s": read_value(
            r"tracing\s+=\s*([0-9.eE+-]+) s", result.stdout, "tracing"
        ),
        "end_to_end_s": read_value(
            r"end-to-end\s+=\s*([0-9.eE+-]+) s", result.stdout, "end-to-end"
        ),
        "field_evaluations": read_integer(
            r"field evaluations =\s*([0-9]+)", result.stdout, "field evaluations"
        ),
        "physical_losses": read_integer(
            r"physical losses =\s*([0-9]+)", result.stdout, "physical losses"
        ),
        "numerical_failures": read_integer(
            r"numerical failures =\s*([0-9]+)",
            result.stdout,
            "numerical failures",
        ),
        "profile": read_profile(profile_path),
    }


def compare(
    candidate: dict[str, object],
    reference: dict[str, object],
    trace_time: float,
) -> dict[str, object]:
    candidate_rows = candidate["profile"]
    reference_rows = reference["profile"]
    assert isinstance(candidate_rows, list)
    assert isinstance(reference_rows, list)
    if len(candidate_rows) != len(reference_rows):
        raise RuntimeError("candidate/reference profile lengths differ")

    state_error = 0.0
    loss_time_error = 0.0
    status_mismatches = 0
    for left, right in zip(candidate_rows, reference_rows):
        assert isinstance(left, dict)
        assert isinstance(right, dict)
        left_lost = float(left["loss_time"]) < 0.999999 * trace_time
        right_lost = float(right["loss_time"]) < 0.999999 * trace_time
        status_mismatches += int(left_lost != right_lost)
        loss_time_error = max(
            loss_time_error,
            abs(float(left["loss_time"]) - float(right["loss_time"])),
        )
        errors = [
            abs(float(left["s_end"]) - float(right["s_end"])),
            angle_difference(float(left["theta_end"]), float(right["theta_end"])),
            angle_difference(float(left["zeta_end"]), float(right["zeta_end"])),
            abs(float(left["speed_end"]) - float(right["speed_end"])),
            abs(float(left["pitch_end"]) - float(right["pitch_end"])),
        ]
        state_error = max(state_error, max(errors))
    return {
        "max_state_error": state_error,
        "max_loss_time_error": loss_time_error,
        "loss_status_mismatches": status_mismatches,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("starts", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--netcdffile", type=Path)
    parser.add_argument("--field-input", type=Path)
    parser.add_argument("--coord-input", type=Path)
    parser.add_argument("--particles", type=int, default=64)
    parser.add_argument("--trace-time", type=float, default=1.0e-4)
    parser.add_argument("--ntimstep", type=int, default=51)
    # This is tight enough to resolve the short production chartmap control
    # without turning the reference itself into the experiment.  A separate
    # convergence spot-check can raise either value when needed.
    parser.add_argument("--reference-npoiper2", type=int, default=256)
    parser.add_argument("--reference-relerr", type=float, default=1.0e-8)
    parser.add_argument("--npoiper2", type=int, nargs="+", default=[32, 64, 128, 256])
    parser.add_argument("--methods", nargs="+", default=["euler", "midpoint"])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.netcdffile is None and (args.field_input is None or args.coord_input is None):
        parser.error("provide --netcdffile or both --field-input and --coord-input")
    if args.work_dir.exists():
        raise RuntimeError(f"work directory already exists: {args.work_dir}")
    args.work_dir.mkdir(parents=True)
    starts_lines = args.starts.read_text().splitlines()
    if len(starts_lines) < args.particles:
        raise RuntimeError("start file has fewer rows than --particles")
    starts = "\n".join(starts_lines[: args.particles]) + "\n"
    executable = args.executable.resolve()
    reference = run_case(
        executable,
        args.work_dir / "reference-dopri",
        starts,
        args.particles,
        "dopri",
        args.reference_npoiper2,
        args.reference_relerr,
        args.trace_time,
        args.ntimstep,
        args.netcdffile,
        args.field_input,
        args.coord_input,
    )
    results: list[dict[str, object]] = []
    for method in args.methods:
        for npoiper2 in args.npoiper2:
            candidate = run_case(
                executable,
                args.work_dir / f"{method}-npoiper{npoiper2}",
                starts,
                args.particles,
                method,
                npoiper2,
                args.reference_relerr,
                args.trace_time,
                args.ntimstep,
                args.netcdffile,
                args.field_input,
                args.coord_input,
            )
            candidate["comparison_to_dopri"] = compare(
                candidate, reference, args.trace_time
            )
            candidate.pop("profile")
            results.append(candidate)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "reference": {key: value for key, value in reference.items() if key != "profile"},
                "results": results,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(json.dumps({"reference": reference["end_to_end_s"], "results": results}, indent=2))


if __name__ == "__main__":
    main()
