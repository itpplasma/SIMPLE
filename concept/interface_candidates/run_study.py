"""Sequential configurations, OpenMP markers; no external SIMPLE dependency."""
from pathlib import Path
import argparse
import gzip
import json
import os
import subprocess
import time

ROOT = Path("concept/interface_candidates")
def configuration(method, case="A", aout=2, speed=.02, tend=4000,
                  tol=1e-9, din=3, order=0, knob=.25, width=0, sample=.5):
    return dict(method=method, case=case, aout=aout, speed=speed, tend=tend,
                tol=tol, din=din, order=order, knob=knob, width=width, sample=sample)

def configurations(pilot=False, horizon=4000, long=False, sampling=False):
    if sampling:
        for case in ("A", "B"):
            base = dict(case=case, aout=2., speed=.02, tend=4000)
            for sample in (.0625, .03125):
                yield configuration("FO", **base, tol=1e-9, sample=sample)
                yield configuration("C3", **base, order=1, sample=sample)
            for tol in (1e-7, 1e-11):
                yield configuration("C3", **base, tol=tol, order=0)
        return
    if long:
        for case, aout in (("A", .75), ("B", 2.)):
            base = dict(case=case, aout=aout, speed=.02, tend=16000)
            for tol in (1e-7, 1e-9, 1e-11):
                yield configuration("FO", **base, tol=tol)
            for order in (0, 1):
                for din in (2, 3, 5, 8):
                    yield configuration("C3", **base, order=order, din=din)
        return
    if pilot:
        for case in ("A", "B"):
            for method, order in (("FO", 0), ("C3", 0), ("C3", 1)):
                yield configuration(method, case=case, order=order, tend=2000)
        return
    for case in ("A", "B"):
        for aout in (.30, .75, 2.):
            for speed in (.01, .02, .04):
                base = dict(case=case, aout=aout, speed=speed, tend=horizon*.02/speed)
                for tol in (1e-7, 1e-9):
                    yield configuration("FO", **base, tol=tol)
                if aout == 2:
                    yield configuration("FO", **base, tol=1e-11)
                for order in (0, 1):
                    for din in (2, 3, 5, 8):
                        yield configuration("C3", **base, order=order, din=din)
                if aout == 2 and speed == .02:
                    for tol in (1e-7, 1e-11):
                        yield configuration("C3", **base, order=1, tol=tol)
                    for sample in (.25, .125):
                        yield configuration("FO", **base, tol=1e-9, sample=sample)
                        yield configuration("C3", **base, order=1, sample=sample)
                    if case == "A":
                        for knob, width in ((.1, .2), (.05, .1), (.025, .05)):
                            yield configuration("GL2", **base, knob=knob, width=width)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--long", action="store_true")
    parser.add_argument("--sampling", action="store_true")
    parser.add_argument("--horizon", type=float, default=4000)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    threads = int(os.environ.get("OMP_NUM_THREADS", "16"))
    if not 1 <= threads <= 16:
        raise SystemExit("OMP_NUM_THREADS must be between 1 and 16")
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ.setdefault("OMP_PROC_BIND", "true")
    os.environ.setdefault("OMP_PLACES", "cores")
    prefix = "pilot" if args.pilot else ("long" if args.long else ("sampling" if args.sampling else "main"))
    results = ROOT / "build" / "study"
    results.mkdir(parents=True, exist_ok=True)
    manifest_path = results / f"{prefix}_manifest.json"
    configurations_list = list(configurations(args.pilot, args.horizon, args.long, args.sampling))
    records = []
    if args.resume and manifest_path.exists():
        records = json.loads(manifest_path.read_text())
    start = time.monotonic()
    for i, cfg in enumerate(configurations_list):
        filename = results / f"{prefix}_{i:03d}.csv"
        if args.resume and i < len(records):
            if records[i]["config"] != cfg:
                raise SystemExit("Resume configuration mismatch; no stale reuse")
            continue
        argv = [str(ROOT / "build" / "compare.x")]
        argv += [str(cfg[k]) for k in ("method", "case", "aout", "speed", "tend",
                                      "tol", "din", "order", "knob", "width", "sample")]
        argv += [str(filename)]
        run = subprocess.run(argv, text=True, capture_output=True, timeout=1800)
        if run.returncode:
            raise SystemExit(f"Failed {cfg}:\n{run.stdout}\n{run.stderr}")
        archive = filename.with_suffix(".csv.gz")
        archive.write_bytes(gzip.compress(filename.read_bytes(), mtime=0))
        filename.unlink()
        records.append(dict(config=cfg, file=str(archive), stdout=run.stdout.strip()))
        manifest_path.write_text(json.dumps(records, indent=2) + "\n")
        elapsed = time.monotonic() - start
        message = f"{prefix} {i+1}/{len(configurations_list)} {run.stdout.strip()} elapsed={elapsed:.1f}s"
        print(message, flush=True)

if __name__ == "__main__":
    main()
