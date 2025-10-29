#!/usr/bin/env python3
"""Parity tests for macrostep orbit output."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pysimple", reason="pysimple module not available")
pytest.importorskip("netCDF4", reason="netCDF4 module not available")

import netCDF4 as nc

import pysimple

REPO_ROOT = Path(__file__).resolve().parents[2]
SIMPLE_EXE = REPO_ROOT / "build" / "simple.x"


def _write_simple_in(path: Path, vmec: Path, n_particles: int, trace_time: float, surface: float) -> None:
    text = f"""&config
startmode = 1
ntestpart = {n_particles}
trace_time = {trace_time:.1e}
npoiper2 = 64
output_orbits_macrostep = .true.
notrace_passing = 0
isw_field_type = 0
integmode = 3
netcdffile = '{vmec}'
deterministic = .true.
sbeg(1) = {surface}
num_surf = 1
/
"""
    (path / "simple.in").write_text(text)


def _read_orbits_nc(path: Path) -> dict[str, np.ndarray]:
    with nc.Dataset(path) as ds:
        return {
            "time": np.array(ds.variables["time"][:]),
            "s": np.array(ds.variables["s"][:]),
            "theta": np.array(ds.variables["theta"][:]),
            "phi": np.array(ds.variables["phi"][:]),
            "p_abs": np.array(ds.variables["p_abs"][:]),
            "v_par": np.array(ds.variables["v_par"][:]),
        }


@pytest.mark.usefixtures("vmec_file")
def test_macrostep_orbit_parity(tmp_path: Path, vmec_file: str) -> None:
    if not SIMPLE_EXE.exists():
        pytest.skip("simple.x is not built; run CMake/Ninja build first")

    vmec_path = Path(vmec_file).resolve()
    n_particles = 4
    surface = 0.35

    fortran_dir = tmp_path / "fortran"
    python_dir = tmp_path / "python"
    fortran_dir.mkdir()
    python_dir.mkdir()

    # Run Fortran executable
    _write_simple_in(fortran_dir, vmec_path, n_particles, 1.0e-3, surface)
    subprocess.run([str(SIMPLE_EXE), "simple.in"], cwd=fortran_dir, check=True)
    assert (fortran_dir / "orbits.nc").exists()

    # Use same start.dat for Python to ensure exact parity
    start_path = python_dir / "start.dat"
    start_path.write_text((fortran_dir / "start.dat").read_text())

    # Run Python API (clean module-level state)
    pysimple.init(
        vmec_path,
        deterministic=True,
        notrace_passing=0,
        npoiper2=64,
        num_surf=1,
        trace_time=1.0e-3,
    )
    # Use wrapper to set sbeg (f90wrap array access workaround)
    pysimple._backend.params_wrapper.set_sbeg(1, surface)

    particles = pysimple.load_particles(start_path)

    # Read Fortran NetCDF output
    fortran_orbits = _read_orbits_nc(fortran_dir / "orbits.nc")

    # Trace particles in Python and collect trajectories
    python_orbits = {
        "time": [],
        "s": [],
        "theta": [],
        "phi": [],
        "p_abs": [],
        "v_par": [],
    }

    for i in range(n_particles):
        result = pysimple.trace_orbit(
            particles[:, i],
            integrator="midpoint",
            return_trajectory=True,
        )

        # Stack trajectories
        python_orbits["time"].append(result["times"])
        python_orbits["s"].append(result["trajectory"][0, :])
        python_orbits["theta"].append(result["trajectory"][1, :])
        python_orbits["phi"].append(result["trajectory"][2, :])
        python_orbits["p_abs"].append(result["trajectory"][3, :])
        python_orbits["v_par"].append(result["trajectory"][4, :])

    # Convert lists to arrays (particle, timestep) to match Fortran output
    for key in python_orbits:
        python_orbits[key] = np.array(python_orbits[key])

    # Compare trajectories between Fortran and Python
    # Note: Python has shape (n_particles, ntimstep), Fortran has (ntimstep, n_particles)
    for key in ["time", "s", "theta", "phi", "p_abs", "v_par"]:
        np.testing.assert_allclose(
            python_orbits[key],
            fortran_orbits[key].T,
            rtol=1e-10,
            atol=1e-12,
            equal_nan=True,  # Allow NaN == NaN for skipped particles
            err_msg=f"Mismatch in {key}",
        )
