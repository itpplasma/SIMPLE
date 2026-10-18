#!/usr/bin/env python3
"""Physical 2 ms cost/loss controls; the phase-consistent matrix is separate."""
import argparse
from pathlib import Path
import subprocess
import sys

parser = argparse.ArgumentParser(__doc__)
parser.add_argument("--output", type=Path, default=Path("build/study/spectre-full"))
parser.add_argument("--markers", type=int, default=256)
parser.add_argument("--duration", type=float, default=2e-3)
parser.add_argument("--threads", type=int, default=16)
args = parser.parse_args()
runner = Path(__file__).with_name("run.py")
for field in ["ncsx", "g3v3"]:
    starts = args.output / f"{field}_fo_default/study_starts.csv"
    for mode, label, grid in [("fo", "default", [32, 64, 32]),
                              ("fo", "fine", [48, 96, 48]),
                              ("hybrid", "default", [32, 64, 32])]:
        output = args.output / f"{field}_{mode}_{label}"
        command = [sys.executable, str(runner), "--output", str(output),
                   "--field", field, "--sampling", "volume", "--mode", mode,
                   "--markers", str(args.markers), "--duration", str(args.duration),
                   "--threads", str(args.threads), "--grid", *map(str, grid)]
        if starts.exists():
            command += ["--starts", str(starts)]
        subprocess.run(command, check=True)
