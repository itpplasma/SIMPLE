#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path


CASES = [
    (
        "cpp6d_vmec_rejected",
        "orbit_model = 5\norbit_coord = 0\n",
        "orbit_model=ORBIT_CPP6D supports only orbit_coord=1",
    ),
    (
        "cp6d_vmec_rejected",
        "orbit_model = 6\norbit_coord = 0\n",
        "orbit_model=ORBIT_CP6D supports only orbit_coord=1",
    ),
    (
        "boris_rejected",
        "orbit_model = 2\n",
        "not available in production alpha-loss tracing",
    ),
    (
        "unknown_rejected",
        "orbit_model = 99\n",
        "unsupported orbit_model",
    ),
]


def write_config(path: Path, body: str) -> None:
    path.write_text("&config\n" + body + "/\n", encoding="utf-8")


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_unsupported_orbit_modes.py SIMPLE_EXE", file=sys.stderr)
        return 2

    simple_exe = Path(sys.argv[1]).resolve()
    workdir = Path.cwd()
    failures = 0

    for name, body, expected in CASES:
        cfg = workdir / f"{name}.in"
        write_config(cfg, body)
        proc = subprocess.run(
            [str(simple_exe), str(cfg)],
            cwd=workdir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        output = proc.stdout + proc.stderr
        if proc.returncode == 0:
            print(f"{name}: expected failure, got success", file=sys.stderr)
            failures += 1
        elif expected not in output:
            print(f"{name}: missing error text: {expected}", file=sys.stderr)
            print(output, file=sys.stderr)
            failures += 1

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
