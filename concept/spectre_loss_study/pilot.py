#!/usr/bin/env python3
"""Sequential paired pilot matrix; completed native cases are resumable."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

parser = argparse.ArgumentParser(__doc__)
parser.add_argument("--baseline", type=Path, required=True)
parser.add_argument("--output", type=Path, default=Path("build/study/spectre"))
parser.add_argument("--duration", type=float, default=2e-6)
parser.add_argument("--markers", type=int, default=16)
parser.add_argument("--threads", type=int, default=8)
parser.add_argument("--variants", nargs="+", choices=["raw", "coarse", "spline", "fine"],
                    default=["raw", "coarse", "spline", "fine"])
args = parser.parse_args()
runner = Path(__file__).with_name("run.py")
for field in ["ncsx", "g3v3"]:
    for sampling in ["sheet", "volume"]:
        starts = args.output / f"{field}_{sampling}_fo_raw/study_starts.csv"
        for mode in ["fo", "hybrid"]:
            for variant, grid in [("raw", None), ("coarse", [16, 32, 16]),
                                  ("spline", [32, 64, 32]),
                                  ("fine", [48, 96, 48])]:
                if variant not in args.variants:
                    continue
                output = args.output / f"{field}_{sampling}_{mode}_{variant}"
                metadata = output / "run.json"
                prior = json.loads(metadata.read_text()) if metadata.exists() else {}
                same = all(prior.get(key) == value for key, value in dict(
                    field=field, sampling=sampling, mode=mode, grid=grid,
                    markers=args.markers, threads=args.threads,
                    duration=args.duration).items())
                root = args.baseline if variant == "raw" else Path.cwd()
                source_paths = ["concept/spectre_loss_study/study_spectre_orbits.f90",
                                "src/spectre_sympl_orbit.f90", "src/spectre_fo_hybrid.f90",
                                "src/spectre_fo_field.f90"]
                hashes = {p: hashlib.sha256((root / p).read_bytes()).hexdigest()
                          for p in source_paths if (root / p).exists()}
                same = same and prior.get("source_hashes") == hashes
                if same and prior.get("returncode") == 0:
                    print(f"reuse {output}", flush=True)
                    continue
                command = [sys.executable, str(runner), "--output", str(output),
                           "--field", field, "--sampling", sampling, "--mode", mode,
                           "--markers", str(args.markers), "--threads", str(args.threads),
                           "--duration", str(args.duration)]
                if variant == "raw":
                    command += ["--root", str(args.baseline)]
                else:
                    command += ["--grid", *map(str, grid)]
                if starts.exists():
                    command += ["--starts", str(starts)]
                subprocess.run(command, check=True)
