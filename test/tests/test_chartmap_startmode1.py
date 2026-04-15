#!/usr/bin/env python3
"""Regression test for chartmap runs with startmode=1."""

from __future__ import annotations

import shutil
import subprocess
import sys
import os
from pathlib import Path

import numpy as np


def write_config(path: Path) -> None:
    path.write_text(
        """&config
trace_time = 1d-6
ntimstep = 4
ntestpart = 8
field_input = 'test_boozer_chartmap.nc'
coord_input = 'test_boozer_chartmap.nc'
isw_field_type = 2
startmode = 1
integmode = 1
relerr = 1d-10
deterministic = .True.
/
""",
        encoding="ascii",
    )


def run_once(workdir: Path, simple_x: Path) -> tuple[np.ndarray, np.ndarray]:
    for name in [
        "start.dat",
        "times_lost.dat",
        "confined_fraction.dat",
        "avg_inverse_t_lost.dat",
        "results.nc",
        "fort.6601",
    ]:
        target = workdir / name
        if target.exists() or target.is_symlink():
            target.unlink()

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    result = subprocess.run(
        [str(simple_x)],
        cwd=workdir,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise SystemExit(f"simple.x failed with exit code {result.returncode}")

    start = np.loadtxt(workdir / "start.dat")
    confined = np.loadtxt(workdir / "confined_fraction.dat")
    if start.ndim == 1:
        start = start[np.newaxis, :]
    if confined.ndim == 1:
        confined = confined[np.newaxis, :]
    return start, confined


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: test_chartmap_startmode1.py <simple.x> <chartmap.nc>")

    simple_x = Path(sys.argv[1]).resolve()
    chartmap = Path(sys.argv[2]).resolve()
    base = Path.cwd() / "chartmap_startmode1_work"

    if base.exists():
        shutil.rmtree(base)
    base.mkdir()
    (base / "test_boozer_chartmap.nc").symlink_to(chartmap)
    write_config(base / "simple.in")
    np.savetxt(
        base / "start.dat",
        np.full((8, 5), 9.0),
        fmt="%.16e",
    )

    start1, confined1 = run_once(base, simple_x)
    start2, confined2 = run_once(base, simple_x)

    assert start1.shape == (8, 5), start1.shape
    assert not np.allclose(start1, 9.0), "chartmap startmode=1 reused an existing start.dat"
    assert np.allclose(start1, start2), "deterministic chartmap starts changed between runs"
    assert np.allclose(confined1, confined2), "deterministic chartmap outputs changed between runs"
    assert np.isclose(confined1[0, 1] + confined1[0, 2], 1.0), confined1[0]
    assert np.count_nonzero(np.abs(start1[:, 1] - start1[0, 1]) > 1.0e-12) > 0, (
        "expected varied sampled poloidal angles"
    )
    assert np.count_nonzero(np.abs(start1[:, 2] - start1[0, 2]) > 1.0e-12) > 0, (
        "expected varied sampled toroidal angles"
    )
    assert np.count_nonzero(np.abs(start1[:, 4] - start1[0, 4]) > 1.0e-12) > 0, (
        "expected varied sampled pitches"
    )


if __name__ == "__main__":
    main()
