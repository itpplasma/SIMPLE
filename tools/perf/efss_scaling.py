#!/usr/bin/env python3
"""Run bounded, reproducible local SIMPLE scaling measurements from the repo root.

Run directories retain the input, stdout and physics outputs. No source changes,
worktrees, downloads, builds or remote execution are performed by this tool.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

import numpy as np


def physics(directory):
    """Read total confinement and unresolved fractions, retaining the time grid."""
    confined = np.loadtxt(directory / "confined_fraction.dat", ndmin=2)
    unresolved = np.loadtxt(directory / "unresolved_fraction.dat", ndmin=2)
    if not np.isfinite(confined).all() or not np.isfinite(unresolved).all():
        raise ValueError(f"Nonfinite fraction output in {directory}")
    if (confined.shape[0] != unresolved.shape[0]
            or not np.array_equal(confined[:, 0], unresolved[:, 0])):
        raise ValueError(f"Fraction time grids differ in {directory}")
    return confined[:, 0], confined[:, 1] + confined[:, 2], unresolved[:, 1]


def check_traced(directory, markers):
    codes = np.loadtxt(directory / "orbit_exit_code.dat", ndmin=2)[:, 1]
    if len(codes) != markers:
        raise ValueError(f"Unexpected marker count in {directory}")
    if np.any(codes == 3):
        raise ValueError(f"Skipped markers invalidate throughput in {directory}")


def phase_seconds(log, gpu):
    pattern = (r"tracing\s*=\s*([\d.]+) s" if gpu else
               r"Particle tracing completed\s+([\d.]+) s")
    match = re.search(pattern, log)
    if match is None:
        raise ValueError("Missing tracing timer in SIMPLE stdout")
    return float(match[1])


def run(args):
    executable, wout = args.executable.resolve(), args.wout.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    metadata = {
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True).strip(),
        "executable": str(executable),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "wout_sha256": hashlib.sha256(wout.read_bytes()).hexdigest(),
        "arguments": {k: str(v) for k, v in vars(args).items()},
    }
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    baseline = None
    rows = []
    # Alternate ascending/descending order to expose time-dependent interference.
    for repetition in range(1, args.repetitions + 1):
        order = args.threads if repetition % 2 else list(reversed(args.threads))
        for threads in order:
            directory = args.output / f"n{args.markers}-t{threads}-r{repetition}"
            directory.mkdir(exist_ok=False)
            (directory / "wout.nc").symlink_to(wout)
            config = f"""&config
 netcdffile = 'wout.nc'
 ntestpart = {args.markers}
 contr_pp = -1d100
 notrace_passing = 0
 sbeg = 0.3d0
 trace_time = {args.trace_time:.12e}
 ntimstep = 1001
 npoiper2 = 128
 isw_field_type = 2
 integmode = 0
 relerr = 1d-9
 deterministic = .true.
 ran_seed = 42
"""
            if args.gpu:
                config += f""" gpu_mode = 'production'
 gpu_backend = 'cuda-native'
 gpu_method = 'dopri'
 gpu_landreman = .true.
 gpu_num_devices = 1
 gpu_t_block = {args.trace_time:.12e}
 gpu_loss_tau = 1d0
 gpu_maxloss = 1d0
 gpu_min_timestep = 0d0
"""
            (directory / "simple.in").write_text(config + "/\n")
            env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_THREAD_LIMIT="16",
                       OMP_DYNAMIC="FALSE", OMP_PLACES="cores", OMP_PROC_BIND="spread",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       CUDA_VISIBLE_DEVICES=str(args.device))
            start = time.perf_counter()
            with (directory / "stdout.log").open("w") as log:
                subprocess.run([str(executable), "simple.in"], cwd=directory,
                               env=env, stdout=log, stderr=subprocess.STDOUT,
                               check=True, timeout=args.timeout)
            wall = time.perf_counter() - start
            log = (directory / "stdout.log").read_text()
            grid, confined, unresolved = physics(directory)
            check_traced(directory, args.markers)
            starts = np.loadtxt(directory / "start.dat", ndmin=2)
            if baseline is None:
                baseline = (grid, confined, unresolved, starts)
            if not np.array_equal(starts, baseline[3]):
                raise ValueError(f"Starting markers changed in {directory}")
            np.testing.assert_allclose(grid, baseline[0], rtol=0, atol=1e-14)
            delta = float(np.max(np.abs(confined - baseline[1])))
            unresolved_delta = float(np.max(np.abs(unresolved - baseline[2])))
            if max(delta, unresolved_delta) > 1e-12:
                raise ValueError(f"Thread-dependent physics in {directory}: {delta}, {unresolved_delta}")
            row = dict(markers=args.markers, threads=threads, repetition=repetition,
                       backend="cuda" if args.gpu else "cpu", wall_s=wall,
                       trace_s=phase_seconds(log, args.gpu),
                       confined_final=float(confined[-1]),
                       unresolved_final=float(unresolved[-1]),
                       max_confined_delta=delta, max_unresolved_delta=unresolved_delta)
            rows.append(row)
            with (args.output / "measurements.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(row))
                writer.writeheader()
                writer.writerows(rows)
            print(json.dumps(row), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=Path("build/simple.x"))
    parser.add_argument("--wout", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markers", type=int, default=256)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8, 16])
    parser.add_argument("--repetitions", type=int, default=2)
    parser.add_argument("--trace-time", type=float, default=0.01)
    parser.add_argument("--timeout", type=float, default=1200)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--device", type=int, default=0)
    args = parser.parse_args()
    if (any(t < 1 or t > 16 for t in args.threads) or args.markers < 1
            or args.repetitions < 1 or args.trace_time <= 0 or args.timeout <= 0):
        parser.error("Use 1..16 threads and positive markers, repetitions and times")
    run(args)


if __name__ == "__main__":
    main()
