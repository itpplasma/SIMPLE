"""Field-level retirement gate for SIMPLE#430: .bc -> boozmn -> chartmap.

Validates that the pre-generated Boozer chartmap (committed as a test fixture)
contains the variables SIMPLE's chartmap reader expects and that Bmod values
agree with direct Fourier summation from the .bc harmonics.

Fixture files committed at test/test_data/:
  circ.bc            -- Strumberger Boozer .bc file (axisymmetric tokamak)
  boozmn_circ.nc     -- booz_xform-format NetCDF generated from circ.bc via
                        libneo bc_to_booz_xform (feat/eqdsk-boozer-chartmap,
                        merged into libneo main as PR #343-345)
  chartmap_circ.nc   -- Boozer chartmap for SIMPLE generated from boozmn_circ.nc
                        via libneo booz_xform_to_boozer_chartmap (nrho=40,
                        ntheta=96, nzeta=1)

No orbit runs; no POTATO; no libneo converters required at test time.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")
_scipy_interp = pytest.importorskip("scipy.interpolate")
CubicSpline = _scipy_interp.CubicSpline
RegularGridInterpolator = _scipy_interp.RegularGridInterpolator

# ---------------------------------------------------------------------------
# Fixture paths: committed into the repo, no external dependencies.
# ---------------------------------------------------------------------------
_TEST_DATA = Path(__file__).resolve().parents[1] / "test_data"
CIRC_BC = _TEST_DATA / "circ.bc"
BOOZMN_NC = _TEST_DATA / "boozmn_circ.nc"
CHARTMAP_NC = _TEST_DATA / "chartmap_circ.nc"

TWOPI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Direct .bc Fourier evaluation (reference)
# ---------------------------------------------------------------------------

def _bmod_from_bc(bc, s_test: np.ndarray, theta_b_test: np.ndarray) -> np.ndarray:
    """Evaluate |B| from .bc Fourier harmonics at (s, theta_B) pairs.

    Axisymmetric: bmnc[m] * cos(m * theta_B).  Spline in s between surfaces.
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
# Tests
# ---------------------------------------------------------------------------

def test_chartmap_has_required_variables():
    """Chartmap file must contain the variables SIMPLE's chartmap reader expects."""
    assert CHARTMAP_NC.exists(), f"committed fixture missing: {CHARTMAP_NC}"
    required = {"rho", "s", "theta", "zeta", "Bmod", "A_phi", "B_theta", "B_phi"}
    with netCDF4.Dataset(CHARTMAP_NC) as ds:
        present = set(ds.variables.keys())
    missing = required - present
    assert not missing, f"chartmap missing variables: {missing}"


def test_bmod_positive():
    """All Bmod values in the chartmap must be positive (field strength)."""
    assert CHARTMAP_NC.exists(), f"committed fixture missing: {CHARTMAP_NC}"
    with netCDF4.Dataset(CHARTMAP_NC) as ds:
        bmod = np.asarray(ds.variables["Bmod"][:], dtype=float)
    assert float(bmod.min()) > 0.0, f"non-positive Bmod found: min={bmod.min():.4g}"


def test_chartmap_bmod_matches_bc_fourier():
    """Bmod in the chartmap agrees with direct .bc Fourier evaluation.

    Conversion chain verified once at fixture-generation time:
        .bc harmonics -> boozmn NetCDF -> splined chartmap (nrho=40, ntheta=96, nzeta=1)

    Expected error budget at mid-radius (0.20 < s < 0.90):
        cubic spline in rho: ~1e-3 to 5e-3 relative
        Fourier truncation: negligible (circ.bc has m0=18 modes)
    Tolerance: 5e-3 relative.
    """
    pytest.importorskip("libneo")
    assert CIRC_BC.exists(), f"committed fixture missing: {CIRC_BC}"
    assert CHARTMAP_NC.exists(), f"committed fixture missing: {CHARTMAP_NC}"

    from libneo.boozer import BoozerFile
    bc = BoozerFile(str(CIRC_BC))

    with netCDF4.Dataset(CHARTMAP_NC) as ds:
        rho_grid = np.asarray(ds.variables["rho"][:], dtype=float)
        theta_grid = np.asarray(ds.variables["theta"][:], dtype=float)
        bmod_nc = np.asarray(ds.variables["Bmod"][:], dtype=float)

    bmod_2d = bmod_nc[0, :, :].T  # (nrho, ntheta)

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

    bmod_direct_T = _bmod_from_bc(bc, s_test, theta_test)

    # chartmap stores Bmod in Gauss; convert to Tesla for comparison
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


def test_bmod_axis_scale():
    """On-axis Bmod (m=0 harmonic) from chartmap matches .bc m=0 mode to 1%."""
    pytest.importorskip("libneo")
    assert CIRC_BC.exists(), f"committed fixture missing: {CIRC_BC}"
    assert CHARTMAP_NC.exists(), f"committed fixture missing: {CHARTMAP_NC}"

    from libneo.boozer import BoozerFile
    bc = BoozerFile(str(CIRC_BC))

    s_surf = np.asarray(bc.s, dtype=float)
    m0 = np.asarray(bc.m[0], dtype=int)
    bmnc = np.array([bc.bmnc[k] for k in range(bc.nsurf)], dtype=float)

    idx_m0 = np.where(m0 == 0)[0]
    assert len(idx_m0) > 0, "no m=0 mode in .bc"
    b00_at_mid = float(CubicSpline(s_surf, bmnc[:, idx_m0[0]])(0.25))

    with netCDF4.Dataset(CHARTMAP_NC) as ds:
        rho_grid = np.asarray(ds.variables["rho"][:], dtype=float)
        bmod_nc = np.asarray(ds.variables["Bmod"][:], dtype=float)

    bmod_2d = bmod_nc[0, :, :].T  # (nrho, ntheta)
    bmod_theta_avg = bmod_2d.mean(axis=1) / 1.0e4  # G -> T

    rho_mid = np.sqrt(0.25)
    b00_chart = float(np.interp(rho_mid, rho_grid, bmod_theta_avg))

    rel_err = abs(b00_chart - b00_at_mid) / abs(b00_at_mid)
    assert rel_err < 1.0e-2, (
        f"on-axis B00 mismatch: bc={b00_at_mid:.6f} T chart={b00_chart:.6f} T "
        f"rel_err={rel_err:.2e}"
    )
