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

import simple

REPO_ROOT = Path(__file__).resolve().parents[2]
SIMPLE_EXE = REPO_ROOT / "build" / "simple.x"


def _write_simple_in(path: Path, vmec: Path, n_particles: int, trace_time: float) -> None:
    text = f"""&config
startmode = 2
ntestpart = {n_particles}
trace_time = {trace_time:.1e}
npoiper2 = 64
output_orbits_macrostep = .true.
notrace_passing = 0
isw_field_type = 0
integmode = 3
netcdffile = '{vmec}'
deterministic = .true.
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
    session = simple.SimpleSession(vmec_path)

    n_particles = 4

    fortran_dir = tmp_path / "fortran"
    python_dir = tmp_path / "python"
    fortran_dir.mkdir()
    python_dir.mkdir()

    _write_simple_in(fortran_dir, vmec_path, n_particles, 1.0e-3)

    with simple.field_type(0):
        with simple.temporary_parameters(deterministic=True):
            batch = session.sample_surface(n_particles, surface=0.35)

        np.savetxt(fortran_dir / "start.dat", batch.positions.T, fmt="%.18e")

        subprocess.run([str(SIMPLE_EXE), "simple.in"], cwd=fortran_dir, check=True)
        assert (fortran_dir / "orbits.nc").exists()

        # Reuse the exact start file consumed by simple.x to drive the Python API
        start_path = python_dir / "start.dat"
        start_path.write_text((fortran_dir / "start.dat").read_text())

        particles = session.load_particles(start_path)
        cwd = os.getcwd()
        try:
            os.chdir(python_dir)
            with simple.temporary_parameters(
                notrace_passing=0,
                npoiper2=64,
                deterministic=True,
            ):
                with simple.macrostep_output(True):
                    results = session.trace(
                        particles, tmax=1.0e-3, integrator="midpoint"
                    )
        finally:
            os.chdir(cwd)

    assert (python_dir / "orbits.nc").exists()
    assert results.n_particles == n_particles

    fortran_data = _read_orbits_nc(fortran_dir / "orbits.nc")
    python_data = _read_orbits_nc(python_dir / "orbits.nc")

    for key in ("time", "s", "theta", "phi", "p_abs", "v_par"):
        mask = ~np.isnan(python_data[key]) & ~np.isnan(fortran_data[key])
        np.testing.assert_allclose(
            python_data[key][mask],
            fortran_data[key][mask],
            rtol=1e-10,
            atol=1e-12,
        )
