#!/usr/bin/env python3
"""Rescale a SPECTRE file and check its field through the Fortran reader."""
import subprocess
import sys
import tempfile
from pathlib import Path

LENGTH, FIELD = 1.7, 18.0


def main(exe, tool, source):
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / "rescaled.h5"
        subprocess.run([sys.executable, tool, source, str(target),
                        "--length", str(LENGTH), "--field", str(FIELD)], check=True)
        subprocess.run([exe, source, str(target), str(LENGTH), str(FIELD)], check=True)


if __name__ == "__main__":
    main(*sys.argv[1:4])
