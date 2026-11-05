#!/usr/bin/env python3
"""Compare the run time of two simple.x builds on short benchmark cases.

Each repetition starts the reference and the candidate at the same moment, one
OpenMP thread each, so both see the same machine load. The decision uses the
per-pair time ratio, separately for preprocessing (everything before the orbit
trace) and for the trace. A phase regresses when the median ratio exceeds
--max-slowdown and every pair is slower.

usage: benchmark.py REF_EXE CUR_EXE [--repeats 5] [--work DIR] [--data DIR]
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import time
import urllib.request

HERE = Path(__file__).resolve().parent
WOUT_URL = ("https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/"
            "wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc")
TRACE = re.compile(r"INFO: Particle tracing completed\s+([0-9.]+) s \(Total:\s+([0-9.]+) s\)")
# Below this reference time the timer and process start dominate.
MIN_PHASE_S = 0.2


def phases(log):
    match = TRACE.search(log.read_text())
    if not match:
        raise RuntimeError(f"no tracing timer in {log}")
    trace, total = map(float, match.groups())
    return {"preprocessing": total - trace, "trace": trace}


def start(exe, case_dir, env):
    log = open(case_dir/"simple.log", "w")
    return subprocess.Popen([str(exe)], cwd=case_dir, env=env, stdout=log,
                            stderr=subprocess.STDOUT), log


def stage(case, run_dir, data):
    run_dir.mkdir(parents=True)
    shutil.copy(HERE/case/"simple.in", run_dir)
    if "netcdffile" in (HERE/case/"simple.in").read_text():
        os.symlink(data/"wout.nc", run_dir/"wout.nc")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("ref", type=Path)
    parser.add_argument("cur", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--max-slowdown", type=float, default=1.2)
    parser.add_argument("--work", type=Path, default=Path("performance-runs"))
    parser.add_argument("--data", type=Path, default=None)
    parser.add_argument("--cases", nargs="*",
                        default=sorted(p.parent.name for p in HERE.glob("*/simple.in")))
    args = parser.parse_args()

    data = args.data or args.work/"data"
    data.mkdir(parents=True, exist_ok=True)
    if not (data/"wout.nc").exists():
        urllib.request.urlretrieve(WOUT_URL, data/"wout.nc.part")
        (data/"wout.nc.part").rename(data/"wout.nc")
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    exes = {"ref": args.ref.resolve(), "cur": args.cur.resolve()}

    times = {}
    for case in args.cases:
        for rep in range(args.repeats):
            dirs = {}
            for build in exes:
                dirs[build] = args.work/case/f"{build}-{rep}"
                shutil.rmtree(dirs[build], ignore_errors=True)
                stage(case, dirs[build], data.resolve())
            started = time.perf_counter()
            procs = {b: start(exes[b], dirs[b], env) for b in exes}
            for build, (proc, log) in procs.items():
                if proc.wait():
                    sys.exit(f"{build} failed in {dirs[build]}/simple.log")
                log.close()
            print(f"{case} repeat {rep}: {time.perf_counter() - started:.1f} s", flush=True)
            for build in exes:
                for phase, value in phases(dirs[build]/"simple.log").items():
                    times.setdefault((case, phase), {}).setdefault(build, []).append(value)

    failed = []
    report = []
    print(f"\n{'case':22s} {'phase':14s} {'ref s':>8s} {'cur s':>8s} {'ratio':>7s}  pair ratios")
    for (case, phase), t in times.items():
        ref, cur = statistics.median(t["ref"]), statistics.median(t["cur"])
        timed = ref >= MIN_PHASE_S
        ratios = [c/r for r, c in zip(t["ref"], t["cur"])] if timed else [float("nan")]
        ratio = statistics.median(ratios)
        slower = timed and ratio > args.max_slowdown and min(ratios) > 1.0
        verdict = "REGRESSION" if slower else ("" if timed else f"not judged (ref < {MIN_PHASE_S} s)")
        print(f"{case:22s} {phase:14s} {ref:8.3f} {cur:8.3f} {ratio:7.3f}  "
              f"{' '.join(f'{x:.3f}' for x in ratios)}  {verdict}")
        report.append(dict(case=case, phase=phase, ref=t["ref"], cur=t["cur"],
                           median_ratio=ratio, regression=slower))
        if slower:
            failed.append(f"{case} {phase}")
    (args.work/"benchmark.json").write_text(json.dumps(dict(
        ref=str(exes["ref"]), cur=str(exes["cur"]), max_slowdown=args.max_slowdown,
        results=report), indent=2) + "\n")
    if failed:
        sys.exit(f"performance regression beyond {args.max_slowdown}x: {', '.join(failed)}")


if __name__ == "__main__":
    main()
