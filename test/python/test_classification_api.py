#!/usr/bin/env python3
"""Tests for the fast classification Python API."""

from __future__ import annotations

import os
import math
import subprocess
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pysimple", reason="pysimple module not available")

import simple
pytest.importorskip("pysimple", reason="pysimple module not available")

REPO_ROOT = Path(__file__).resolve().parents[2]
SIMPLE_EXE = REPO_ROOT / "build" / "simple.x"


@pytest.mark.usefixtures("vmec_file")
def test_classify_fast_restores_parameters(vmec_file: str) -> None:
    session = simple.SimpleSession(vmec_file)
    batch = session.sample_surface(6, surface=0.35)

    params_before = simple.get_parameters("tcut", "fast_class", "class_plot", "trace_time")

    result = session.classify_fast(batch)

    params_after = simple.get_parameters("tcut", "fast_class", "class_plot", "trace_time")
    assert params_before == params_after

    assert result.j_parallel.shape == (batch.n_particles,)
    assert result.topology.shape == (batch.n_particles,)
    assert result.minkowski is None


def _load_fort_table(path: Path) -> np.ndarray:
    if not path.exists() or path.stat().st_size == 0:
        return np.empty((0, 3))
    data = np.loadtxt(path)
    return np.atleast_2d(data)


def test_classify_fast_matches_fortran(tmp_path: Path, vmec_file: str) -> None:
    if not SIMPLE_EXE.exists():
        pytest.skip("simple.x is not built; run CMake/Ninja build first")

    vmec_path = Path(vmec_file).resolve()
    session = simple.SimpleSession(vmec_path)

    n_particles = 8
    batch = session.sample_surface(n_particles, surface=0.4)

    fortran_dir = tmp_path / "fortran"
    python_dir = tmp_path / "python"
    fortran_dir.mkdir()
    python_dir.mkdir()

    np.savetxt(fortran_dir / "start.dat", batch.positions.T, fmt="%.18e")
    np.savetxt(python_dir / "start.dat", batch.positions.T, fmt="%.18e")

    config_text = f"""&config
startmode = 1
ntestpart = {n_particles}
trace_time = 1.0d-1
fast_class = .true.
class_plot = .true.
tcut = 1.0d-1
notrace_passing = 1
isw_field_type = 2
deterministic = .true.
netcdffile = '{vmec_path}'
/
"""
    config_path = fortran_dir / "simple.in"
    config_path.write_text(config_text)

    subprocess.run([str(SIMPLE_EXE), config_path.name], cwd=fortran_dir, check=True)

    particles = session.load_particles(python_dir / "start.dat")
    original_params = simple.get_parameters("notrace_passing", "deterministic", "trace_time")
    cwd = os.getcwd()
    try:
        os.chdir(python_dir)
        simple.set_parameters(notrace_passing=1, deterministic=True, trace_time=0.1)
        result = session.classify_fast(
            particles,
            classification_time=0.1,
            legacy_files=True,
            output_dir=Path.cwd(),
        )
        start_snapshot = simple._backend.snapshot_start_positions(n_particles)
    finally:
        os.chdir(cwd)
        simple.set_parameters(**original_params)
    theta = np.mod(start_snapshot[1, :], 2 * math.pi)
    pitch = start_snapshot[4, :]
    trap = result.trap_parameter

    def match_row(row: np.ndarray) -> int:
        theta_row = row[0] % (2 * math.pi)
        angle_diff = np.abs(((theta - theta_row + math.pi) % (2 * math.pi)) - math.pi)
        total = angle_diff + np.abs(pitch - row[1]) + np.abs(trap - row[2])
        idx = int(np.argmin(total))
        assert total[idx] < 1e-5, "Unable to match classification record"
        return idx

    matched = set()

    def check_unit(unit: int, code: int, target: np.ndarray) -> None:
        python_rows = _load_fort_table(python_dir / f"fort.{unit}")
        if python_rows.size == 0:
            return
        fortran_rows = _load_fort_table(fortran_dir / f"fort.{unit}")
        assert fortran_rows.size > 0, f"Fortran output missing for fort.{unit}"
        remaining = fortran_rows.copy()
        for row in python_rows:
            diff = np.max(np.abs(remaining - row), axis=1)
            best = int(np.argmin(diff))
            assert diff[best] < 1e-6, f"Fortran output missing row from fort.{unit}"
            idx = match_row(row)
            matched.add(idx)
            assert target[idx] == code, (
                f"Mismatch for particle {idx}: expected code {code}, got {target[idx]}"
            )
            remaining = np.delete(remaining, best, axis=0)

    check_unit(40012, 1, result.j_parallel)
    check_unit(40022, 2, result.j_parallel)
    check_unit(40032, 0, result.j_parallel)

    check_unit(50012, 1, result.topology)
    check_unit(50022, 2, result.topology)
    check_unit(50032, 0, result.topology)
