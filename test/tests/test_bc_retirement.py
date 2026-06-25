#!/usr/bin/env python3
"""Field-level retirement gate for SIMPLE#430: .bc -> boozmn -> chartmap.

Generates a reference Boozer chartmap from circ.bc using the libneo converters
  bc_to_booz_xform  (feat/eqdsk-boozer-chartmap worktree)
  booz_xform_to_boozer_chartmap  (same worktree)

then checks that the Bmod values stored in the chartmap (which is exactly what
SIMPLE reads at runtime via its chartmap I/O module) match direct Fourier
summation from the .bc harmonics to within the spline interpolation error.
No orbit runs; no POTATO; no .bc reader in SIMPLE (SIMPLE has none).

Usage::
    pytest test_bc_retirement.py -v
    ctest -R test_bc_retirement
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import netCDF4
import numpy as np
import pytest
from scipy.interpolate import CubicSpline, RegularGridInterpolator

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CIRC_BC = Path("/home/ert/code/NEO-RT/examples/circ.bc")
EQDSK_BOOZ_DIR = Path("/tmp/booz-wt/eqdsk-booz/python/libneo")
TWOPI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Load converters from the feat/eqdsk-boozer-chartmap worktree.
# The installed libneo package at /home/ert/code/libneo/python does not yet
# contain bc_to_booz_xform; load those modules directly by file path so
# we do not depend on a specific installation order.
# ---------------------------------------------------------------------------

def _load_converter(name: str):
    """Return a module loaded from EQDSK_BOOZ_DIR/<name>.py."""
    fullname = f"libneo.{name}"
    if fullname in sys.modules:
        return sys.modules[fullname]
    path = EQDSK_BOOZ_DIR / f"{name}.py"
    spec = importlib.util.spec_from_file_location(fullname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[fullname] = mod
    spec.loader.exec_module(mod)
    return mod


def _ensure_converters():
    for name in ("boozer_chartmap_writer", "booz_xform_to_boozer_chartmap",
                 "bc_to_booz_xform"):
        _load_converter(name)


# ---------------------------------------------------------------------------
# Direct .bc Fourier evaluation (reference)
# ---------------------------------------------------------------------------

def _bmod_from_bc(bc, s_test: np.ndarray, theta_b_test: np.ndarray) -> np.ndarray:
    """Evaluate |B| directly from .bc Fourier harmonics at (s, theta_B) pairs.

    For axisymmetric equilibria (n0b=0): bmnc[m] cos(m theta_B).
    Uses cubic spline in s between .bc surfaces.
    """
    nsurf = bc.nsurf
    s_surf = np.asarray(bc.s, dtype=float)
    m0 = np.asarray(bc.m[0], dtype=int)
    bmnc = np.array([bc.bmnc[k] for k in range(nsurf)], dtype=float)

    n_test = len(s_test)
    bmod = np.empty(n_test, dtype=float)
    for i in range(n_test):
        bmnc_at_s = np.array(
            [float(CubicSpline(s_surf, bmnc[:, j])(s_test[i])) for j in range(len(m0))],
            dtype=float,
        )
        bmod[i] = float(np.dot(bmnc_at_s, np.cos(m0 * theta_b_test[i])))
    return bmod


# ---------------------------------------------------------------------------
# Pytest fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def circ_bc_file():
    if not CIRC_BC.exists():
        pytest.skip(f".bc file not found: {CIRC_BC}")
    pytest.importorskip("netCDF4")
    pytest.importorskip("scipy")
    if not EQDSK_BOOZ_DIR.exists():
        pytest.skip(f"eqdsk-booz worktree not found: {EQDSK_BOOZ_DIR}")
    return CIRC_BC


@pytest.fixture(scope="module")
def chartmap_path(tmp_path_factory, circ_bc_file):
    """Build circ.bc -> boozmn -> chartmap once per test session."""
    _ensure_converters()
    bc_to_booz_xform = sys.modules["libneo.bc_to_booz_xform"]
    booz_xform_to_boozer_chartmap = sys.modules["libneo.booz_xform_to_boozer_chartmap"]

    tmp = tmp_path_factory.mktemp("bc_retirement")
    boozmn = tmp / "boozmn_circ.nc"
    chartmap = tmp / "chartmap_circ.nc"

    bc_to_booz_xform.convert_bc_to_boozmn(circ_bc_file, boozmn)
    assert boozmn.exists(), "bc_to_booz_xform produced no output"

    booz_xform_to_boozer_chartmap.convert_boozmn_to_chartmap(
        boozmn, chartmap, nrho=40, ntheta=96, nzeta=1
    )
    assert chartmap.exists(), "booz_xform_to_boozer_chartmap produced no output"
    return chartmap


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_chartmap_has_required_variables(chartmap_path):
    """Chartmap file must contain the variables SIMPLE's chartmap reader expects."""
    required = {"rho", "s", "theta", "zeta", "Bmod", "A_phi", "B_theta", "B_phi"}
    with netCDF4.Dataset(chartmap_path) as ds:
        present = set(ds.variables.keys())
    missing = required - present
    assert not missing, f"chartmap missing variables: {missing}"


def test_bmod_positive(chartmap_path):
    """All Bmod values in the chartmap must be positive (field strength)."""
    with netCDF4.Dataset(chartmap_path) as ds:
        bmod = np.asarray(ds.variables["Bmod"][:], dtype=float)
    assert float(bmod.min()) > 0.0, f"non-positive Bmod found: min={bmod.min():.4g}"


def test_chartmap_bmod_matches_bc_fourier(chartmap_path, circ_bc_file):
    """Bmod in the chartmap agrees with direct .bc Fourier evaluation.

    The conversion chain is:
        .bc harmonics -> boozmn NetCDF -> splined chartmap (nrho=40, ntheta=96, nzeta=1)

    Expected error budget at mid-radius (0.20 < s < 0.90):
        cubic spline in rho: ~1e-3 to 5e-3 relative
        Fourier truncation: negligible (circ.bc has m0=18 modes)
    Tolerance: 5e-3 relative (conservative, matches libneo's own cross-path test).
    """
    from libneo.boozer import BoozerFile

    bc = BoozerFile(str(circ_bc_file))

    with netCDF4.Dataset(chartmap_path) as ds:
        rho_grid = np.asarray(ds.variables["rho"][:], dtype=float)
        theta_grid = np.asarray(ds.variables["theta"][:], dtype=float)
        # Bmod shape: (nzeta, ntheta, nrho); pick zeta=0 slice
        bmod_nc = np.asarray(ds.variables["Bmod"][:], dtype=float)

    bmod_2d = bmod_nc[0, :, :].T  # shape (nrho, ntheta)

    # Periodic wrap for theta interpolation.
    theta_w = np.append(theta_grid, TWOPI)
    bmod_w = np.concatenate([bmod_2d, bmod_2d[:, :1]], axis=1)
    interp = RegularGridInterpolator(
        (rho_grid, theta_w),
        bmod_w,
        method="cubic",
        bounds_error=False,
        fill_value=None,
    )

    rng = np.random.default_rng(42)
    n_test = 30
    s_test = rng.uniform(0.20, 0.90, n_test)
    rho_test = np.sqrt(s_test)
    theta_test = rng.uniform(0.0, TWOPI, n_test)

    # Direct Fourier evaluation in Tesla from the .bc file.
    bmod_direct_T = _bmod_from_bc(bc, s_test, theta_test)

    # Chartmap stores Bmod in Gauss; convert to Tesla for comparison.
    bmod_chart_G = interp(np.column_stack([rho_test, theta_test]))
    bmod_chart_T = bmod_chart_G / 1.0e4

    tol = 5.0e-3
    for i in range(n_test):
        rel_err = abs(bmod_direct_T[i] - bmod_chart_T[i]) / abs(bmod_direct_T[i])
        assert rel_err < tol, (
            f"point {i}: s={s_test[i]:.3f} theta={theta_test[i]:.4f} rad "
            f"Bmod_bc={bmod_direct_T[i]:.6f} T "
            f"Bmod_chart={bmod_chart_T[i]:.6f} T "
            f"rel_err={rel_err:.2e} > tol={tol:.0e}"
        )


def test_bmod_axis_scale(chartmap_path, circ_bc_file):
    """On-axis Bmod (m=0, n=0 harmonic) from chartmap matches .bc m=0 mode to 1 %."""
    from libneo.boozer import BoozerFile

    bc = BoozerFile(str(circ_bc_file))
    s_surf = np.asarray(bc.s, dtype=float)
    m0 = np.asarray(bc.m[0], dtype=int)
    bmnc = np.array([bc.bmnc[k] for k in range(bc.nsurf)], dtype=float)
    # m=0 mode coefficient at mid-radius (s=0.25)
    idx_m0 = np.where(m0 == 0)[0]
    assert len(idx_m0) > 0, "no m=0 mode in .bc"
    b00_at_mid = float(CubicSpline(s_surf, bmnc[:, idx_m0[0]])(0.25))

    with netCDF4.Dataset(chartmap_path) as ds:
        rho_grid = np.asarray(ds.variables["rho"][:], dtype=float)
        theta_grid = np.asarray(ds.variables["theta"][:], dtype=float)
        bmod_nc = np.asarray(ds.variables["Bmod"][:], dtype=float)

    bmod_2d = bmod_nc[0, :, :].T  # (nrho, ntheta)
    # Average over theta to get m=0 content.
    bmod_theta_avg = bmod_2d.mean(axis=1) / 1.0e4  # convert G -> T

    rho_mid = np.sqrt(0.25)
    b00_chart = float(np.interp(rho_mid, rho_grid, bmod_theta_avg))

    rel_err = abs(b00_chart - b00_at_mid) / abs(b00_at_mid)
    assert rel_err < 1.0e-2, (
        f"on-axis B00 mismatch: bc={b00_at_mid:.6f} T chart={b00_chart:.6f} T "
        f"rel_err={rel_err:.2e}"
    )
