#!/usr/bin/env python3
"""Check native segmented loss accounting against hand-constructed outcomes."""

from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def extract(pattern: str, output: str) -> str:
    match = re.search(pattern, output)
    if match is None:
        raise AssertionError(f"missing output matching {pattern!r}")
    return match.group(1)


def run_case(
    executable: Path,
    wout: Path,
    starts_source: Path,
    method: str,
    relerr: float,
) -> None:
    with tempfile.TemporaryDirectory(prefix=f"landreman-{method}-", dir=Path.cwd()) as tmp:
        root = Path(tmp)
        starts = []
        boundary_count = 8 if method == "dopri" else 4
        for index, line in enumerate(starts_source.read_text().splitlines()):
            columns = line.split()
            if index < boundary_count:
                columns[0] = "1.0"
            starts.append(" ".join(columns))
        (root / "start.dat").write_text("\n".join(starts) + "\n")
        (root / "simple.in").write_text(
            "&config\n"
            "trace_time = 1d-4\n"
            "ntimstep = 101\n"
            "ntestpart = 8\n"
            "startmode = 2\n"
            f"netcdffile = '{wout}'\n"
            "isw_field_type = 2\n"
            "integmode = 1\n"
            "npoiper2 = 256\n"
            f"relerr = {relerr:.17e}\n"
            "/\n"
        )
        environment = os.environ.copy()
        environment.update(
            {
                "OMP_NUM_THREADS": "1",
                "SIMPLE_GPU_RUN": "1",
                "SIMPLE_GPU_METHOD": method,
                "SIMPLE_GPU_LANDREMAN": "1",
                "SIMPLE_GPU_T_BLOCK": "1e-6",
                "SIMPLE_GPU_MAXLOSS": "0.1",
                "SIMPLE_GPU_LOSS_TAU": "0.1",
                "SIMPLE_GPU_MIN_TIMESTEP": "0",
            }
        )
        result = subprocess.run(
            [executable, "simple.in"],
            cwd=root,
            env=environment,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        output = result.stdout
        losses = int(extract(r"physical losses =\s*([0-9]+)", output))
        failures = int(extract(r"numerical failures =\s*([0-9]+)", output))
        observed = float(extract(r"observed until =\s*([0-9.eE+-]+)", output))
        energy_loss = float(
            extract(r"energy loss fraction =\s*([0-9.eE+-]+)", output)
        )
        if losses != boundary_count or failures != 0:
            raise AssertionError(
                f"{method}: expected {boundary_count} physical losses and no failures, "
                f"got {losses} and {failures}"
            )
        if abs(observed - 1.0e-6) > 5.0e-12:
            raise AssertionError(f"{method}: early stop occurred at {observed}")
        expected_loss = boundary_count / 8.0
        if abs(energy_loss - expected_loss) > 0.01:
            raise AssertionError(f"{method}: energy-loss fraction is {energy_loss}")

        times = []
        for line in (root / "times_lost.dat").read_text().splitlines():
            columns = line.split()
            times.append(float(columns[1]))
        if len(times) != 8:
            raise AssertionError(f"{method}: expected 8 loss times, got {len(times)}")
        if not all(0.0 <= value < 1.0e-6 for value in times[:boundary_count]):
            raise AssertionError(
                f"{method}: boundary particles were not lost: "
                f"{times[:boundary_count]}"
            )
        if not all(value == 1.0e-4 for value in times[boundary_count:]):
            raise AssertionError(
                f"{method}: early-stop survivors were not reported at tmax: "
                f"{times[boundary_count:]}"
            )


def run_short_dopri_trace(
    executable: Path, wout: Path, starts_source: Path
) -> None:
    """Exercise the composed RK and Boozer device path, not only accounting."""
    with tempfile.TemporaryDirectory(prefix="landreman-dopri-orbit-", dir=Path.cwd()) as tmp:
        root = Path(tmp)
        (root / "start.dat").write_text(starts_source.read_text())
        (root / "simple.in").write_text(
            "&config\n"
            "trace_time = 1d-10\n"
            "ntimstep = 2\n"
            "ntestpart = 8\n"
            "startmode = 2\n"
            f"netcdffile = '{wout}'\n"
            "isw_field_type = 2\n"
            "integmode = 1\n"
            "npoiper2 = 256\n"
            "relerr = 1d-6\n"
            "/\n"
        )
        environment = os.environ.copy()
        environment.update(
            {
                "OMP_NUM_THREADS": "1",
                "SIMPLE_GPU_RUN": "1",
                "SIMPLE_GPU_METHOD": "dopri",
                "SIMPLE_GPU_LANDREMAN": "1",
                "SIMPLE_GPU_T_BLOCK": "1e-10",
                "SIMPLE_GPU_MAXLOSS": "1",
                "SIMPLE_GPU_LOSS_TAU": "0.1",
                "SIMPLE_GPU_MIN_TIMESTEP": "1e-12",
            }
        )
        output = subprocess.run(
            [executable, "simple.in"],
            cwd=root,
            env=environment,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        ).stdout
        losses = int(extract(r"physical losses =\s*([0-9]+)", output))
        failures = int(extract(r"numerical failures =\s*([0-9]+)", output))
        evaluations = int(extract(r"field evaluations =\s*([0-9]+)", output))
        if losses != 0 or failures != 0:
            raise AssertionError(
                "short DOPRI orbit must remain confined without numerical failures"
            )
        if evaluations != 288:
            raise AssertionError(
                f"short DOPRI orbit used {evaluations} rather than 288 field evaluations"
            )


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: test_gpu_landreman_segments.py SIMPLE WOUT STARTS")
    executable, wout, starts = (Path(value).resolve() for value in sys.argv[1:])
    run_short_dopri_trace(executable, wout, starts)
    run_case(executable, wout, starts, "dopri", 1.0e-6)
    run_case(executable, wout, starts, "euler", 1.0e-13)
    print("segmented Landreman RK and symplectic accounting passed")


if __name__ == "__main__":
    main()
