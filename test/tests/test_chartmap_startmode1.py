#!/usr/bin/env python3
"""Regression test for chartmap runs with startmode=1."""

from __future__ import annotations

import shutil
import subprocess
import sys
import os
from pathlib import Path

import numpy as np


def write_config(
    path: Path,
    *,
    grid_density: float | None = None,
    generate_start_only: bool = False,
) -> None:
    grid_line = "" if grid_density is None else f"grid_density = {grid_density:.16e}\n"
    start_only_line = (
        "generate_start_only = .True.\n" if generate_start_only else ""
    )
    path.write_text(
        f"""&config
trace_time = 1d-6
ntimstep = 4
ntestpart = 8
sbeg = 0.5d0
{grid_line}field_input = 'test_boozer_chartmap.nc'
coord_input = 'test_boozer_chartmap.nc'
isw_field_type = 2
startmode = 1
integmode = 1
relerr = 1d-10
deterministic = .True.
{start_only_line}/
""",
        encoding="ascii",
    )


def run_once(
    workdir: Path,
    simple_x: Path,
    *,
    expect_orbit_output: bool = True,
) -> tuple[np.ndarray, np.ndarray | None]:
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
    if start.ndim == 1:
        start = start[np.newaxis, :]

    confined_path = workdir / "confined_fraction.dat"
    if not expect_orbit_output:
        assert not confined_path.exists(), "generate_start_only traced chartmap orbits"
        assert not (workdir / "times_lost.dat").exists(), (
            "generate_start_only wrote chartmap orbit output"
        )
        return start, None

    confined = np.loadtxt(confined_path)
    if confined.ndim == 1:
        confined = confined[np.newaxis, :]
    return start, confined


def run_case(
    workdir: Path,
    simple_x: Path,
    *,
    grid_density: float | None = None,
    generate_start_only: bool = False,
) -> tuple[np.ndarray, np.ndarray | None]:
    write_config(
        workdir / "simple.in",
        grid_density=grid_density,
        generate_start_only=generate_start_only,
    )
    np.savetxt(
        workdir / "start.dat",
        np.full((8, 5), 9.0),
        fmt="%.16e",
    )
    start, confined = run_once(
        workdir,
        simple_x,
        expect_orbit_output=not generate_start_only,
    )
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

    start1, confined1 = run_case(base, simple_x)
    start2, confined2 = run_case(base, simple_x)

    assert confined1 is not None
    assert confined2 is not None
    assert start1.shape == (8, 5), start1.shape
    assert not np.allclose(start1, 9.0), "chartmap startmode=1 reused an existing start.dat"
    np.testing.assert_allclose(start1[:, 0], np.sqrt(0.5), rtol=0.0, atol=1.0e-15)
    assert np.array_equal(start1, start2), (
        "deterministic chartmap starts changed between runs"
    )
    assert np.array_equal(confined1, confined2), (
        "deterministic chartmap outputs changed between runs"
    )
    assert np.isclose(confined1[0, 1] + confined1[0, 2], 1.0, rtol=0.0, atol=1.0e-15), confined1[0]
    assert np.count_nonzero(np.abs(start1[:, 1] - start1[0, 1]) > 1.0e-12) > 0, (
        "expected varied sampled poloidal angles"
    )
    assert np.count_nonzero(np.abs(start1[:, 2] - start1[0, 2]) > 1.0e-12) > 0, (
        "expected varied sampled toroidal angles"
    )
    assert np.count_nonzero(np.abs(start1[:, 4] - start1[0, 4]) > 1.0e-12) > 0, (
        "expected varied sampled pitches"
    )

    grid_start, _ = run_case(base, simple_x, grid_density=0.25)
    assert grid_start.shape == (9, 5), grid_start.shape
    np.testing.assert_allclose(grid_start[:, 0], np.sqrt(0.5), rtol=0.0, atol=1.0e-15)

    start_only, confined = run_case(base, simple_x, generate_start_only=True)
    assert confined is None
    np.testing.assert_allclose(start_only[:, 0], np.sqrt(0.5), rtol=0.0, atol=1.0e-15)


if __name__ == "__main__":
    main()
