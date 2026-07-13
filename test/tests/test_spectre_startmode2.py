#!/usr/bin/env python3
"""SPECTRE startmode=2 must preserve explicit starting points."""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


STARTS = np.array([
    [0.25, 0.1, 0.2, 1.0, 0.3],
    [1.25, 1.1, 0.7, 2.0, -0.2],
    [2.25, 2.2, 0.4, 0.5, 0.8],
])


def config(field):
    return f"""&config
field_input = '{field}'
integ_coords = 6
integmode = 0
ntestpart = {len(STARTS)}
startmode = 2
generate_start_only = .True.
spectre_ncon_r = 8
spectre_ncon_th = 8
spectre_ncon_phi = 8
/
"""


def main():
    binary, field = map(Path, sys.argv[1:])
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        np.savetxt(work / "start.dat", STARTS, fmt="%.17g")
        (work / "simple.in").write_text(config(field.resolve()))
        proc = subprocess.run([binary.resolve(), "simple.in"], cwd=work,
                              capture_output=True, text=True, timeout=120)
        if proc.returncode != 0:
            raise RuntimeError(proc.stdout + proc.stderr)
        actual = np.loadtxt(work / "start.dat")
        if not np.array_equal(actual, STARTS):
            raise AssertionError(f"explicit starts were replaced:\n{actual}")
    print("SPECTRE startmode=2 explicit starts PASS")


if __name__ == "__main__":
    main()
