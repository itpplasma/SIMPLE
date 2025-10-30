#!/usr/bin/env python3
"""Parity tests for classification between Python and Fortran."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pysimple", reason="pysimple module not available")

import pysimple

REPO_ROOT = Path(__file__).resolve().parents[2]
SIMPLE_EXE = REPO_ROOT / "build" / "simple.x"


def _write_simple_in(path: Path, vmec: Path, n_particles: int, trace_time: float, surface: float) -> None:
    text = f"""&config
startmode = 1
ntestpart = {n_particles}
trace_time = {trace_time:.1e}
npoiper2 = 64
tcut = 0.1
class_plot = .false.
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


def _read_fortran_classification(path: Path, n_particles: int) -> dict[str, np.ndarray]:
    """Read classification results from Fortran times_lost.dat."""
    data = np.loadtxt(path / "times_lost.dat")

    # times_lost.dat format:
    # ipart, passing(0/1), loss_time, trap_par, perp_inv
    # Classification is in iclass array but not written to file!
    # We'll rely on passing/trapped and loss times for basic parity

    return {
        "passing": data[:, 1].astype(bool),
        "loss_times": data[:, 2],
        "trap_parameter": data[:, 3],
        "perpendicular_invariant": data[:, 4],
    }


@pytest.mark.usefixtures("vmec_file")
def test_classification_parity(tmp_path: Path, vmec_file: str) -> None:
    """Compare classification results between Fortran and Python."""
    if not SIMPLE_EXE.exists():
        pytest.skip("simple.x is not built; run CMake/Ninja build first")

    vmec_path = Path(vmec_file).resolve()
    n_particles = 10
    surface = 0.35

    fortran_dir = tmp_path / "fortran"
    python_dir = tmp_path / "python"
    fortran_dir.mkdir()
    python_dir.mkdir()

    # Run Fortran executable
    _write_simple_in(fortran_dir, vmec_path, n_particles, 1.0e-3, surface)
    subprocess.run([str(SIMPLE_EXE), "simple.in"], cwd=fortran_dir, check=True)
    assert (fortran_dir / "times_lost.dat").exists()

    # Use same start.dat for Python to ensure exact parity
    start_path = python_dir / "start.dat"
    start_path.write_text((fortran_dir / "start.dat").read_text())

    # Run Python API
    pysimple.init(
        vmec_path,
        deterministic=True,
        notrace_passing=0,
        npoiper2=64,
        num_surf=1,
        trace_time=1.0e-3,
        tcut=0.1,
        class_plot=False,
    )
    # Note: sbeg is set by sampling, not manually

    particles = pysimple.load_particles(start_path)
    python_results = pysimple.classify_parallel(particles, integrator="midpoint")

    # Read Fortran results
    fortran_results = _read_fortran_classification(fortran_dir, n_particles)

    # Compare basic results (classification indices not written to file)
    np.testing.assert_array_equal(
        python_results["passing"],
        fortran_results["passing"],
        err_msg="Passing classification mismatch",
    )

    np.testing.assert_allclose(
        python_results["loss_times"],
        fortran_results["loss_times"],
        rtol=1e-10,
        atol=1e-12,
        equal_nan=True,
        err_msg="Loss times mismatch",
    )

    np.testing.assert_allclose(
        python_results["trap_parameter"],
        fortran_results["trap_parameter"],
        rtol=1e-10,
        atol=1e-12,
        equal_nan=True,
        err_msg="Trap parameter mismatch",
    )

    np.testing.assert_allclose(
        python_results["perpendicular_invariant"],
        fortran_results["perpendicular_invariant"],
        rtol=1e-10,
        atol=1e-12,
        equal_nan=True,
        err_msg="Perpendicular invariant mismatch",
    )

    # Verify classification arrays are populated (values depend on ntcut/class_plot settings)
    assert "minkowski" in python_results
    assert "jpar" in python_results
    assert "topology" in python_results
    assert python_results["minkowski"].shape == (n_particles,)
    assert python_results["jpar"].shape == (n_particles,)
    assert python_results["topology"].shape == (n_particles,)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
