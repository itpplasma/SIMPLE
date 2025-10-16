#!/usr/bin/env python3
"""Tests for the fast classification Python API."""

from __future__ import annotations

import os
import math
import subprocess
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pytest

import simple

REPO_ROOT = Path(__file__).resolve().parents[2]
SIMPLE_EXE = REPO_ROOT / "build" / "simple.x"


@pytest.mark.usefixtures("vmec_file")
def test_classify_fast_restores_parameters(vmec_file: str) -> None:
    simple.load_vmec(vmec_file)
    batch = simple.ParticleBatch(6)
    batch.initialize_from_samplers(vmec_file, method="surface", s=0.35)

    params_before = simple.get_parameters("tcut", "fast_class", "class_plot", "trace_time")

    result = simple.classify_fast(batch)

    params_after = simple.get_parameters("tcut", "fast_class", "class_plot", "trace_time")
    assert params_before == params_after

    assert result.j_parallel.shape == (batch.n_particles,)
    assert result.topology.shape == (batch.n_particles,)
    assert result.minkowski.shape == (batch.n_particles,)


def _load_fort_table(path: Path) -> np.ndarray:
    if not path.exists() or path.stat().st_size == 0:
        return np.empty((0, 3))
    data = np.loadtxt(path)
    return np.atleast_2d(data)


def test_classify_fast_matches_fortran(tmp_path: Path, vmec_file: str) -> None:
    if not SIMPLE_EXE.exists():
        pytest.skip("simple.x is not built; run CMake/Ninja build first")

    vmec_path = Path(vmec_file).resolve()
    simple.load_vmec(vmec_path)

    n_particles = 8
    batch = simple.ParticleBatch(n_particles)
    batch.initialize_from_samplers(vmec_path, method="surface", s=0.4)

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

    particles = simple.load_particle_file(vmec_path, python_dir / "start.dat")
    original_params = simple.get_parameters("notrace_passing", "deterministic", "trace_time")
    cwd = os.getcwd()
    try:
        os.chdir(python_dir)
        simple.load_vmec(vmec_path, force=True)
        simple.set_parameters(notrace_passing=1, deterministic=True, trace_time=0.1)
        result = simple.classify_fast(particles, tcut=0.1, vmec_file=vmec_path, write_files=True)
        start_snapshot = simple._backend.snapshot_start_positions(n_particles)
    finally:
        os.chdir(cwd)
        simple.set_parameters(**original_params)
    def fort_count(unit: int) -> int:
        return _load_fort_table(fortran_dir / f"fort.{unit}").shape[0]

    j_counts = {
        0: fort_count(40032),
        1: fort_count(40012),
        2: fort_count(40022),
    }
    topo_counts = {
        0: fort_count(50032),
        1: fort_count(50012),
        2: fort_count(50022),
    }

    for code, expected in j_counts.items():
        actual = int(np.count_nonzero(result.j_parallel == code))
        if expected > 0:
            assert actual > 0, f"Python classification missing J_parallel code {code}"

    for code, expected in topo_counts.items():
        actual = int(np.count_nonzero(result.topology == code))
        if expected > 0:
            assert actual > 0, f"Python classification missing topology code {code}"
