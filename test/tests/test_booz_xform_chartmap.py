#!/usr/bin/env python3
"""Validate the booz_xform -> Boozer chartmap converter on the QA equilibrium.

Three layers, all against independent references:
  1. Spectral reference with grid convergence: cubic interpolation of the
     tabulated chartmap must approach direct summation of the boozmn
     harmonics as the chartmap grid is refined.
  2. Golden cross-check: the converted chartmap must agree with the chartmap
     written by SIMPLE's own export_boozer_chartmap from the same wout.
  3. End-to-end orbits: confined fractions of a SIMPLE run on the converted
     chartmap must match the VMEC-direct run within the 1/npart resolution,
     and the derived rmajor must keep dtaumin within 5 percent.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import netCDF4
import numpy as np
from scipy.interpolate import CubicSpline, RegularGridInterpolator

from boozer_chartmap_artifacts import (
    confined_metrics,
    download_if_missing,
    load_confined_fraction,
    run_cmd,
)

SCRIPT_DIR = Path(__file__).resolve().parent
BUILD_DIR = SCRIPT_DIR.parent.parent / "build"
SIMPLE_X = BUILD_DIR / "simple.x"
TOOL_X = BUILD_DIR / "test" / "tests" / "export_boozer_chartmap_tool.x"
TOOLS_DIR = SCRIPT_DIR.parent.parent / "tools"
CONVERTER = TOOLS_DIR / "booz_xform_to_boozer_chartmap.py"
TEST_DATA = SCRIPT_DIR.parent / "test_data"
WOUT = TEST_DATA / "wout.nc"
WOUT_URL = (
    "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/"
    "tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
)

sys.path.insert(0, str(TOOLS_DIR))
from booz_xform_to_boozer_chartmap import (  # noqa: E402
    fourier_eval,
    interp_coeffs,
    read_boozmn,
)

TWOPI = 2.0 * np.pi


def fail(msg: str) -> None:
    raise SystemExit(f"FAIL: {msg}")


def check(label: str, value: float, tol: float) -> None:
    status = "ok" if value <= tol else "FAIL"
    print(f"  {label}: {value:.3e} (tol {tol:.0e}) {status}")
    if value > tol:
        fail(f"{label} = {value:.3e} exceeds {tol:.0e}")


def make_boozmn(workdir: Path) -> Path:
    import booz_xform

    boozmn = workdir / "boozmn_qa.nc"
    b = booz_xform.Booz_xform()
    b.read_wout(str(WOUT))
    b.mboz = 48
    b.nboz = 48
    b.verbose = 0
    b.run()
    b.write_boozmn(str(boozmn))
    return boozmn


def convert(boozmn: Path, output: Path, nrho: int, ntheta: int, nzeta: int) -> Path:
    run_cmd(
        [
            sys.executable,
            str(CONVERTER),
            str(boozmn),
            str(output),
            "--nrho", str(nrho),
            "--ntheta", str(ntheta),
            "--nzeta", str(nzeta),
        ],
        cwd=output.parent,
        label=f"convert {output.name}",
    )
    return output


def load_chartmap(path: Path) -> dict:
    with netCDF4.Dataset(path) as ds:
        out = {
            "rho": ds.variables["rho"][:].data,
            "s": ds.variables["s"][:].data,
            "theta": ds.variables["theta"][:].data,
            "zeta": ds.variables["zeta"][:].data,
            "A_phi": ds.variables["A_phi"][:].data,
            "B_theta": ds.variables["B_theta"][:].data,
            "B_phi": ds.variables["B_phi"][:].data,
            "Bmod": np.transpose(ds.variables["Bmod"][:].data, (2, 1, 0)),
            "x": np.transpose(ds.variables["x"][:].data, (2, 1, 0)),
            "y": np.transpose(ds.variables["y"][:].data, (2, 1, 0)),
            "z": np.transpose(ds.variables["z"][:].data, (2, 1, 0)),
            "nfp": int(ds.variables["num_field_periods"][:]),
            "torflux": float(ds.torflux),
            "attrs": list(ds.ncattrs()),
        }
    return out


def bmod_interpolator(c: dict) -> RegularGridInterpolator:
    bmod = c["Bmod"]
    nrho, n_th, n_ze = bmod.shape
    bmod_periodic = np.empty((nrho, n_th + 1, n_ze + 1))
    bmod_periodic[:, :n_th, :n_ze] = bmod
    bmod_periodic[:, n_th, :n_ze] = bmod[:, 0, :]
    bmod_periodic[:, :n_th, n_ze] = bmod[:, :, 0]
    bmod_periodic[:, n_th, n_ze] = bmod[:, 0, 0]
    theta_f = np.linspace(0.0, TWOPI, n_th + 1)
    zeta_f = np.linspace(0.0, TWOPI / c["nfp"], n_ze + 1)
    return RegularGridInterpolator((c["rho"], theta_f, zeta_f), bmod_periodic,
                                   method="cubic")


def spectral_bmod(d: dict, pts: np.ndarray) -> np.ndarray:
    """Direct boozmn summation at (rho, theta, zeta) points, in Gauss."""
    rho_half = np.sqrt(d["s_half"])
    out = np.empty(len(pts))
    for i, (rho, th, ze) in enumerate(pts):
        c = interp_coeffs(d["bmnc"], d["ixm"], rho_half, np.array([rho]))
        val = fourier_eval(c, d["ixm"], d["ixn"], np.array([th]),
                           np.array([ze]), "cos")[0, 0, 0]
        out[i] = val * 1.0e4
    return out


def check_spectral_convergence(d: dict, workdir: Path, boozmn: Path) -> None:
    print("== Spectral reference with grid convergence ==")
    coarse = convert(boozmn, workdir / "chartmap_coarse.nc", 25, 16, 32)
    fine = convert(boozmn, workdir / "chartmap_fine.nc", 50, 48, 96)

    rng = np.random.default_rng(20260610)
    n_pts = 200
    pts = np.column_stack([
        rng.uniform(0.2, 0.9, n_pts),
        rng.uniform(0.0, TWOPI, n_pts),
        rng.uniform(0.0, TWOPI / d["nfp"], n_pts),
    ])
    ref = spectral_bmod(d, pts)

    errs = {}
    for name, path in (("coarse", coarse), ("fine", fine)):
        c = load_chartmap(path)
        approx = bmod_interpolator(c)(pts)
        errs[name] = float(np.max(np.abs(approx - ref) / np.abs(ref)))
        print(f"  Bmod {name} grid max rel err vs spectral: {errs[name]:.3e}")

    if errs["fine"] >= 0.5 * errs["coarse"]:
        fail(
            f"no grid convergence: fine err {errs['fine']:.3e} "
            f">= 0.5 * coarse err {errs['coarse']:.3e}"
        )
    check("Bmod fine-grid spectral error", errs["fine"], 1.0e-4)


def check_golden_export(workdir: Path, chartmap: Path) -> Path:
    print("== Golden cross-check against export_boozer_chartmap ==")
    export_dir = workdir / "export"
    export_dir.mkdir(exist_ok=True)
    start_vmec = export_dir / "start_vmec.dat"
    start_vmec.write_text("0.3 0.5 0.5 1.0 0.1\n", encoding="utf-8")
    ref_chartmap = export_dir / "chartmap_vmec.nc"
    run_cmd(
        [
            str(TOOL_X),
            str(WOUT),
            str(ref_chartmap),
            str(start_vmec),
            str(export_dir / "start_boozer.dat"),
        ],
        cwd=export_dir,
        label="export_boozer_chartmap",
    )

    ref = load_chartmap(ref_chartmap)
    new = load_chartmap(chartmap)
    for c, name in ((ref, "exported"), (new, "converted")):
        if "rmajor" in c["attrs"]:
            fail(f"{name} chartmap still carries the removed rmajor attribute")

    check("torflux rel diff",
          abs(new["torflux"] - ref["torflux"]) / abs(ref["torflux"]), 1.0e-6)

    b_scale = float(np.max(np.abs(ref["B_phi"])))
    interp = CubicSpline(new["s"], new["A_phi"])(ref["s"])
    check("A_phi max diff / scale",
          float(np.max(np.abs(interp - ref["A_phi"]))
                / float(np.max(np.abs(ref["A_phi"])))), 1.0e-5)

    for q, scale, tol in (
        ("B_theta", b_scale, 1.0e-4),
        ("B_phi", b_scale, 1.0e-4),
    ):
        interp = CubicSpline(new["rho"], new[q])(ref["rho"])
        check(f"{q} max diff / scale",
              float(np.max(np.abs(interp - ref[q])) / scale), tol)

    rng = np.random.default_rng(42)
    n_pts = 400
    pts = np.column_stack([
        rng.uniform(0.15, 0.95, n_pts),
        rng.uniform(0.0, TWOPI, n_pts),
        rng.uniform(0.0, TWOPI / ref["nfp"], n_pts),
    ])
    b_ref = bmod_interpolator(ref)(pts)
    b_new = bmod_interpolator(new)(pts)
    check("Bmod max rel diff",
          float(np.max(np.abs(b_new - b_ref) / np.abs(b_ref))), 5.0e-3)

    geom_new = {
        q: RegularGridInterpolator(
            (new["rho"], new["theta"], new["zeta"]), new[q],
            method="cubic", bounds_error=False, fill_value=None)
        for q in ("x", "y", "z")
    }
    th_grid, ze_grid = np.meshgrid(ref["theta"], ref["zeta"], indexing="ij")
    max_diff = 0.0
    for ir in (5, 25, 45):
        p = np.column_stack([
            np.full(th_grid.size, ref["rho"][ir]),
            th_grid.ravel(),
            ze_grid.ravel(),
        ])
        for q in ("x", "y", "z"):
            diff = np.max(np.abs(geom_new[q](p) - ref[q][ir].ravel()))
            max_diff = max(max_diff, float(diff))
    check("geometry max |diff| [cm]", max_diff, 1.0)

    # Derived major radius (axis-average R) must agree between the two
    # writers; it intentionally differs from the wout volume-based Rmajor_p.
    r_ref = float(np.mean(np.sqrt(ref["x"][0]**2 + ref["y"][0]**2))) / 100.0
    r_new = float(np.mean(np.sqrt(new["x"][0]**2 + new["y"][0]**2))) / 100.0
    check("derived rmajor rel diff", abs(r_new - r_ref) / r_ref, 1.0e-3)
    return export_dir / "start_boozer.dat"


SIMPLE_IN_VMEC = """\
&config
multharm = 5
contr_pp = -1e10
trace_time = 1d-2
macrostep_time_grid = 'log'
ntimstep = 61
sbeg = 0.3d0
ntestpart = 32
netcdffile = 'wout.nc'
isw_field_type = 2
deterministic = .True.
integmode = 1
facE_al = 1.0d0
/
"""

SIMPLE_IN_CHARTMAP = """\
&config
multharm = 5
contr_pp = -1e10
trace_time = 1d-2
macrostep_time_grid = 'log'
ntimstep = 61
sbeg = 0.3d0
ntestpart = 32
field_input = 'boozer_chartmap.nc'
coord_input = 'boozer_chartmap.nc'
isw_field_type = 2
deterministic = .True.
integmode = 1
startmode = 2
facE_al = 1.0d0
/
"""


def parse_tau_line(stdout: str) -> tuple[float, int]:
    for line in stdout.splitlines():
        if line.strip().startswith("tau:"):
            parts = line.split()
            return float(parts[2].replace("D", "E")), int(parts[4])
    fail("run did not print a tau line")


def check_e2e_orbits(workdir: Path, chartmap: Path) -> None:
    print("== End-to-end VMEC vs booz_xform chartmap orbits ==")
    vmec_dir = workdir / "vmec_run"
    chart_dir = workdir / "chartmap_run"
    for run_dir in (vmec_dir, chart_dir):
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir()

    shutil.copy2(WOUT, vmec_dir / "wout.nc")
    (vmec_dir / "simple.in").write_text(SIMPLE_IN_VMEC, encoding="utf-8")
    vmec_stdout = run_cmd([str(SIMPLE_X)], cwd=vmec_dir, label="VMEC run")

    # Convert the VMEC start.dat to Boozer/chartmap reference coordinates.
    start_boozer = chart_dir / "start.dat"
    run_cmd(
        [
            str(TOOL_X),
            str(WOUT),
            str(chart_dir / "chartmap_unused.nc"),
            str(vmec_dir / "start.dat"),
            str(start_boozer),
        ],
        cwd=chart_dir,
        label="start.dat conversion",
    )
    (chart_dir / "chartmap_unused.nc").unlink()
    shutil.copy2(chartmap, chart_dir / "boozer_chartmap.nc")
    (chart_dir / "simple.in").write_text(SIMPLE_IN_CHARTMAP, encoding="utf-8")
    chart_stdout = run_cmd([str(SIMPLE_X)], cwd=chart_dir, label="chartmap run")

    # Microstep: rmajor is now the geometry-derived axis-average R, which for
    # QA differs from the wout volume-based Rmajor_p by 1.2 percent; dtaumin
    # and ntau inherit that difference linearly. 5 percent bounds it with
    # margin while still catching order-of-magnitude rmajor bugs (350b0f2).
    tau_tol = 5.0e-2
    ref_dtaumin, ref_ntau = parse_tau_line(vmec_stdout)
    new_dtaumin, new_ntau = parse_tau_line(chart_stdout)
    check("dtaumin rel diff",
          abs(new_dtaumin - ref_dtaumin) / abs(ref_dtaumin), tau_tol)
    check("ntau rel diff", abs(new_ntau - ref_ntau) / ref_ntau, tau_tol)

    ref = confined_metrics(load_confined_fraction(vmec_dir / "confined_fraction.dat"))
    new = confined_metrics(load_confined_fraction(chart_dir / "confined_fraction.dat"))
    tol = 1.0 / ref["npart"]
    for key in ("total_conf", "pass_conf", "trap_conf"):
        n = min(len(ref[key]), len(new[key]))
        diff = float(np.max(np.abs(np.asarray(ref[key][:n]) - np.asarray(new[key][:n]))))
        check(f"confined fraction max |{key} diff|", diff, tol)


def main() -> None:
    for path, label in ((SIMPLE_X, "simple.x"), (TOOL_X, "export tool"),
                        (CONVERTER, "converter script")):
        if not path.exists():
            raise SystemExit(f"Missing {label}: {path}")
    download_if_missing(WOUT, WOUT_URL)

    workdir = Path.cwd() / "booz_xform_chartmap_test"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)

    print("== Running booz_xform on QA wout ==")
    boozmn = make_boozmn(workdir)
    d = read_boozmn(boozmn)

    check_spectral_convergence(d, workdir, boozmn)
    chartmap = workdir / "chartmap_fine.nc"
    check_golden_export(workdir, chartmap)
    check_e2e_orbits(workdir, chartmap)
    print("BOOZ_XFORM CHARTMAP CONVERTER TEST PASSED")


if __name__ == "__main__":
    main()
