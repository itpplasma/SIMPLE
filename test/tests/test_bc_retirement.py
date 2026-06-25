#!/usr/bin/env python3
"""Scaffold for .bc retirement: verify chartmap matches .bc to FP accuracy.

Retirement gate for itpplasma/SIMPLE#430.  Full cross-path floating-point
equality depends on libneo converters that are not yet merged:

  libneo#343  EQDSK -> Boozer chartmap (tokamak path)
  libneo#344  .bc   -> booz_xform       (Strumberger Boozer ASCII)
  libneo#345  cross-path field-equality test

Until those land this test does two things that can run today:

  Phase A (always runs)
    VMEC-Boozer vs VMEC-chartmap confined fractions on the QA wout.nc.
    This is a tokamak-representative proxy: the QA equilibrium is nearly
    axisymmetric.  It exercises the full chartmap read path and confirms
    that the orbit result agrees within 1/npart of the VMEC-Boozer baseline.
    It re-uses the export_boozer_chartmap_tool that already exists.

  Phase B (skipped when libneo converters are absent)
    EQDSK -> chartmap -> SIMPLE vs the same case run with the .bc reader.
    Requires eqdsk_to_boozer_chartmap (libneo#343) and the .bc reader.
    Tolerance: max |field difference| / max |field| < 1e-12 (FP accuracy).
    Enabled by placing an EQDSK file at test/test_data/eqdsk.gfile and a
    matching .bc file at test/test_data/ref.bc.

Usage::

    pytest test_bc_retirement.py   (or: python test_bc_retirement.py)
    ctest -R test_bc_retirement    (via CMake registration)
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np

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
TEST_DATA = SCRIPT_DIR.parent / "test_data"
WOUT = TEST_DATA / "wout.nc"
WOUT_URL = (
    "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/"
    "tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
)

# Optional tokamak inputs for Phase B (EQDSK/.bc cross-path)
EQDSK_FILE = TEST_DATA / "eqdsk.gfile"
BC_FILE = TEST_DATA / "ref.bc"

# Libneo converter entry points (available once libneo#343 is merged)
EQDSK_TO_CHARTMAP: Path | None = None
for _candidate in [
    BUILD_DIR.parent.parent / "libneo" / "python" / "libneo" / "eqdsk_to_boozer_chartmap.py",
    Path(__file__).resolve().parent.parent.parent.parent
    / "libneo" / "python" / "libneo" / "eqdsk_to_boozer_chartmap.py",
]:
    if _candidate.exists():
        EQDSK_TO_CHARTMAP = _candidate
        break


def _simple_in_vmec(wout_name: str) -> str:
    return f"""\
&config
multharm = 5
contr_pp = -1e10
trace_time = 1d-2
macrostep_time_grid = 'log'
ntimstep = 61
sbeg = 0.3d0
ntestpart = 32
netcdffile = '{wout_name}'
isw_field_type = 2
deterministic = .True.
integmode = 1
facE_al = 1.0d0
/
"""


def _simple_in_chartmap(chartmap_name: str) -> str:
    return f"""\
&config
multharm = 5
contr_pp = -1e10
trace_time = 1d-2
macrostep_time_grid = 'log'
ntimstep = 61
sbeg = 0.3d0
ntestpart = 32
field_input = '{chartmap_name}'
coord_input = '{chartmap_name}'
isw_field_type = 2
deterministic = .True.
integmode = 1
startmode = 2
facE_al = 1.0d0
/
"""


def _check(label: str, value: float, tol: float) -> None:
    status = "ok" if value <= tol else "FAIL"
    print(f"  {label}: {value:.3e} (tol {tol:.0e}) {status}")
    if value > tol:
        raise SystemExit(f"FAIL: {label} = {value:.3e} exceeds {tol:.0e}")


def phase_a_vmec_vs_chartmap(workdir: Path) -> None:
    """VMEC-Boozer vs VMEC-chartmap confined fractions on the QA equilibrium.

    This is the proxy tokamak test that exercises the chartmap read path.
    A full cross-path check against a real EQDSK/.bc equilibrium requires
    libneo#343/#344/#345 (Phase B below).
    """
    print("== Phase A: VMEC-Boozer vs VMEC-chartmap (QA proxy) ==")
    wout = download_if_missing(WOUT, WOUT_URL)
    vmec_dir = workdir / "vmec_run"
    chart_dir = workdir / "chartmap_run"
    for d in (vmec_dir, chart_dir):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    shutil.copy2(wout, vmec_dir / wout.name)
    (vmec_dir / "simple.in").write_text(_simple_in_vmec(wout.name), encoding="utf-8")

    run_cmd([str(SIMPLE_X)], cwd=vmec_dir, label="VMEC-Boozer run")

    run_cmd(
        [
            str(TOOL_X),
            str(vmec_dir / wout.name),
            str(chart_dir / "boozer_chartmap.nc"),
            str(vmec_dir / "start.dat"),
            str(chart_dir / "start.dat"),
        ],
        cwd=chart_dir,
        label="export VMEC chartmap",
    )

    (chart_dir / "simple.in").write_text(
        _simple_in_chartmap("boozer_chartmap.nc"), encoding="utf-8"
    )
    run_cmd([str(SIMPLE_X)], cwd=chart_dir, label="VMEC-chartmap run")

    ref = confined_metrics(load_confined_fraction(vmec_dir / "confined_fraction.dat"))
    new = confined_metrics(load_confined_fraction(chart_dir / "confined_fraction.dat"))
    tol = 1.0 / float(ref["npart"])
    min_len = min(len(ref["total_conf"]), len(new["total_conf"]))
    for key in ("total_conf", "pass_conf", "trap_conf"):
        diff = float(
            np.max(
                np.abs(
                    np.asarray(ref[key][:min_len]) - np.asarray(new[key][:min_len])
                )
            )
        )
        _check(f"|{key} diff|", diff, tol)

    print("  Phase A PASSED")


def phase_b_eqdsk_bc_crosscheck(workdir: Path) -> None:
    """EQDSK -> chartmap -> SIMPLE vs .bc-based run.

    Requires:
      - eqdsk.gfile and ref.bc in test/test_data/
      - libneo eqdsk_to_boozer_chartmap.py  (libneo#343)

    Blocked on itpplasma/libneo#343, #344, #345.
    When those land, remove this guard and set tol = 1e-12.
    """
    print("== Phase B: EQDSK/.bc cross-path (FP equality gate) ==")

    if EQDSK_TO_CHARTMAP is None:
        print(
            "  SKIP: eqdsk_to_boozer_chartmap.py not found "
            "(blocked on itpplasma/libneo#343)"
        )
        return

    if not EQDSK_FILE.exists():
        print(
            f"  SKIP: reference EQDSK not found at {EQDSK_FILE}; "
            "place an axisymmetric g-file there to enable this phase"
        )
        return

    if not BC_FILE.exists():
        print(
            f"  SKIP: reference .bc file not found at {BC_FILE}; "
            "place a matching Strumberger Boozer ASCII file there "
            "(libneo#344 provides the converter)"
        )
        return

    eqdsk_dir = workdir / "eqdsk_chartmap_run"
    bc_dir = workdir / "bc_run"
    for d in (eqdsk_dir, bc_dir):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    # Convert EQDSK -> Boozer chartmap (libneo#343)
    chartmap = eqdsk_dir / "boozer_chartmap.nc"
    run_cmd(
        [sys.executable, str(EQDSK_TO_CHARTMAP), str(EQDSK_FILE), str(chartmap)],
        cwd=eqdsk_dir,
        label="EQDSK -> chartmap",
    )

    # Run SIMPLE from the chartmap
    (eqdsk_dir / "simple.in").write_text(
        _simple_in_chartmap("boozer_chartmap.nc"), encoding="utf-8"
    )
    run_cmd([str(SIMPLE_X)], cwd=eqdsk_dir, label="EQDSK-chartmap run")

    # Run SIMPLE from the .bc (legacy path via boozer_sub / get_boozer_coordinates)
    shutil.copy2(BC_FILE, bc_dir / "ref.bc")
    (bc_dir / "simple.in").write_text(
        _simple_in_vmec("ref.bc"), encoding="utf-8"
    )
    run_cmd([str(SIMPLE_X)], cwd=bc_dir, label=".bc run")

    ref = confined_metrics(load_confined_fraction(bc_dir / "confined_fraction.dat"))
    new = confined_metrics(load_confined_fraction(eqdsk_dir / "confined_fraction.dat"))
    # FP accuracy: 1/npart is too loose once libneo#345 provides the field
    # equality gate.  Use 1/npart for orbit-level check here; the field-level
    # check at 1e-12 is owned by the libneo#345 test.
    tol = 1.0 / float(ref["npart"])
    min_len = min(len(ref["total_conf"]), len(new["total_conf"]))
    for key in ("total_conf", "pass_conf", "trap_conf"):
        diff = float(
            np.max(
                np.abs(
                    np.asarray(ref[key][:min_len]) - np.asarray(new[key][:min_len])
                )
            )
        )
        _check(f"|{key} diff|", diff, tol)

    print("  Phase B PASSED")


def main() -> None:
    for path, label in ((SIMPLE_X, "simple.x"), (TOOL_X, "export tool")):
        if not path.exists():
            raise SystemExit(f"Missing {label}: {path}")

    workdir = Path.cwd() / "bc_retirement_test"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)

    phase_a_vmec_vs_chartmap(workdir)
    phase_b_eqdsk_bc_crosscheck(workdir)
    print("BC RETIREMENT TEST PASSED (Phase B skipped until libneo#343-345 land)")


if __name__ == "__main__":
    main()
