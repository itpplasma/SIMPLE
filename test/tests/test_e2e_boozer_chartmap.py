#!/usr/bin/env python3
"""End-to-end benchmark for VMEC and GVEC Boozer chartmaps on common cases."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import netCDF4
import numpy as np

from boozer_chartmap_artifacts import (
    confined_metrics,
    download_if_missing,
    load_confined_fraction,
    plot_all_cases_total_losses,
    plot_case_loss_comparison,
    plot_field_component_comparison,
    plot_surface_comparison,
    summarize_confined_fraction,
)


SCRIPT_DIR = Path(__file__).resolve().parent
BUILD_DIR = SCRIPT_DIR.parent.parent / "build"
SIMPLE_X = BUILD_DIR / "simple.x"
TOOL_X = BUILD_DIR / "test" / "tests" / "export_boozer_chartmap_tool.x"
ROUNDTRIP_X = BUILD_DIR / "test" / "tests" / "test_boozer_chartmap_roundtrip.x"
ROUNDTRIP_PLOT = SCRIPT_DIR / "plot_boozer_chartmap_roundtrip.py"
GVEC_TO_CHARTMAP = SCRIPT_DIR / "gvec_to_boozer_chartmap.py"
TEST_DATA = SCRIPT_DIR.parent / "test_data"


EQUILIBRIA = [
    {
        "name": "QA",
        "wout": TEST_DATA / "wout.nc",
        "url": "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc",
        "trace_time": "1d-2",
        "facE_al": "1.0d0",
        "use_gvec": False,
        "vmec_roundtrip_field_tol": "1.0e-4",
        "gvec_roundtrip_field_tol": "2.5e-4",
        "gvec_minimize_tol": 1.0e-10,
        "gvec_max_iter": 4,
        "gvec_total_iter": 8,
    },
    {
        "name": "QH",
        "wout": TEST_DATA / "wout_qh.nc",
        "url": "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc",
        "trace_time": "1d-2",
        "facE_al": "1.0d0",
        "use_gvec": False,
        "vmec_roundtrip_field_tol": "1.0e-4",
        "gvec_minimize_tol": 1.0e-10,
        "gvec_max_iter": 12,
        "gvec_total_iter": 24,
    },
    {
        "name": "NCSX",
        "wout": TEST_DATA / "wout_ncsx.nc",
        "url": "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/tests/test_files/wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc",
        "trace_time": "1d-2",
        "facE_al": "100.0d0",
        "use_gvec": False,
        "vmec_roundtrip_field_tol": "1.0e-4",
        "gvec_minimize_tol": 1.0e-10,
        "gvec_max_iter": 8,
        "gvec_total_iter": 16,
    },
]


def simple_in_vmec(wout_name: str, trace_time: str, facE_al: str) -> str:
    return f"""\
&config
multharm = 5
contr_pp = -1e10
trace_time = {trace_time}
macrostep_time_grid = 'log'
ntimstep = 61
sbeg = 0.3d0
ntestpart = 32
netcdffile = '{wout_name}'
isw_field_type = 2
deterministic = .True.
integmode = 1
facE_al = {facE_al}
/
"""


def simple_in_chartmap(chartmap_name: str, trace_time: str, facE_al: str) -> str:
    return f"""\
&config
multharm = 5
contr_pp = -1e10
trace_time = {trace_time}
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
facE_al = {facE_al}
/
"""


def run_cmd(cmd: list[str], cwd: Path, label: str, timeout: int = 3600) -> None:
    print(f"[{label}] {' '.join(cmd)}")
    result = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    if result.returncode == 0:
        return
    print(result.stdout[-4000:])
    print(result.stderr[-4000:])
    raise RuntimeError(f"{label} failed with exit code {result.returncode}")


def build_gvec_chartmap(
    wout: Path,
    run_dir: Path,
    minimize_tol: float,
    max_iter: int,
    total_iter: int,
) -> Path:
    import gvec
    from gvec.scripts.convert_wout import convert_vmec_wout
    from gvec.util import read_parameters

    convert_dir = run_dir / "convert"
    solve_dir = run_dir / "solve"
    convert_dir.mkdir(parents=True, exist_ok=True)
    solve_dir.mkdir(parents=True, exist_ok=True)
    convert_vmec_wout(wout, convert_dir)

    restart_state = next(convert_dir.glob("*State*.dat"))
    params = read_parameters(convert_dir / "parameter.toml")
    params["minimize_tol"] = float(minimize_tol)
    params["maxIter"] = int(max_iter)
    params["totalIter"] = int(total_iter)
    gvec.run(
        params,
        restartstate=restart_state,
        runpath=solve_dir,
        redirect_gvec_stdout=True,
        quiet=True,
    )

    param_final = sorted(solve_dir.glob("parameter*_final.ini"))[-1]
    state_final = sorted(solve_dir.glob("*State_final.dat"))[-1]
    output = run_dir / "boozer_chartmap.nc"
    run_cmd(
        [
            sys.executable,
            str(GVEC_TO_CHARTMAP),
            str(param_final),
            str(state_final),
            str(output),
            "--nrho",
            "50",
            "--ntheta",
            "36",
            "--nphi",
            "81",
        ],
        cwd=run_dir,
        label=f"{wout.stem} GVEC->chartmap",
        timeout=3600,
    )
    return output


def run_roundtrip(
    wout: Path,
    chartmap: Path,
    mode: str,
    prefix: Path,
    title: str,
    chart_label: str,
    field_tol: str,
    orbit_tol: str,
) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    run_cmd(
        [
            str(ROUNDTRIP_X),
            str(wout),
            str(chartmap),
            mode,
            str(prefix),
            field_tol,
            orbit_tol,
        ],
        cwd=prefix.parent,
        label=title,
        timeout=3600,
    )
    run_cmd(
        [
            sys.executable,
            str(ROUNDTRIP_PLOT),
            "--prefix",
            str(prefix),
            "--direct-label",
            "VMEC-Boozer",
            "--chartmap-label",
            chart_label,
            "--title",
            title,
        ],
        cwd=prefix.parent,
        label=f"{title} plots",
    )


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def assert_boozer_chartmap_file(path: Path) -> None:
    with netCDF4.Dataset(path) as ds:
        if int(getattr(ds, "boozer_field")) != 1:
            raise RuntimeError(f"{path}: boozer_field attribute missing or invalid")
        if str(getattr(ds, "zeta_convention")) != "boozer":
            raise RuntimeError(f"{path}: expected zeta_convention='boozer'")
        if str(getattr(ds, "rho_convention")) != "rho_tor":
            raise RuntimeError(f"{path}: expected rho_convention='rho_tor'")


def run_single_equilibrium(eq: dict[str, str | Path], outdir: Path) -> tuple[bool, dict[str, dict[str, float]]]:
    name = str(eq["name"])
    wout = download_if_missing(Path(eq["wout"]), str(eq["url"]))
    trace_time = str(eq["trace_time"])
    facE_al = str(eq["facE_al"])
    use_gvec = bool(eq.get("use_gvec", False))
    vmec_roundtrip_field_tol = str(eq.get("vmec_roundtrip_field_tol", "1.0e-4"))
    gvec_roundtrip_field_tol = str(eq.get("gvec_roundtrip_field_tol", "2.5e-4"))
    gvec_minimize_tol = float(eq.get("gvec_minimize_tol", 1.0e-10))
    gvec_max_iter = int(eq.get("gvec_max_iter", 4))
    gvec_total_iter = int(eq.get("gvec_total_iter", 8))

    eq_dir = outdir / name
    if eq_dir.exists():
        shutil.rmtree(eq_dir)
    vmec_dir = eq_dir / "vmec_run"
    chartmap_dir = eq_dir / "vmec_chartmap_run"
    gvec_dir = eq_dir / "gvec_chartmap_run"
    vmec_dir.mkdir(parents=True)
    chartmap_dir.mkdir()
    gvec_dir.mkdir()

    wout_name = wout.name
    shutil.copy2(wout, vmec_dir / wout_name)
    (vmec_dir / "simple.in").write_text(simple_in_vmec(wout_name, trace_time, facE_al), encoding="utf-8")

    print(f"[{name}] Step 1: VMEC Boozer run")
    run_cmd([str(SIMPLE_X)], cwd=vmec_dir, label=f"{name} VMEC")
    cf_vmec = load_confined_fraction(vmec_dir / "confined_fraction.dat")
    start_vmec = vmec_dir / "start.dat"

    print(f"[{name}] Step 2: export VMEC Boozer chartmap")
    chartmap_nc = chartmap_dir / "boozer_chartmap.nc"
    start_boozer = chartmap_dir / "start.dat"
    run_cmd(
        [
            str(TOOL_X),
            str(vmec_dir / wout_name),
            str(chartmap_nc),
            str(start_vmec),
            str(start_boozer),
        ],
        cwd=chartmap_dir,
        label=f"{name} export",
    )
    assert_boozer_chartmap_file(chartmap_nc)

    print(f"[{name}] Step 3: VMEC chartmap SIMPLE run")
    (chartmap_dir / "simple.in").write_text(
        simple_in_chartmap("boozer_chartmap.nc", trace_time, facE_al),
        encoding="utf-8",
    )
    run_cmd([str(SIMPLE_X)], cwd=chartmap_dir, label=f"{name} VMEC chartmap")
    cf_chartmap = load_confined_fraction(chartmap_dir / "confined_fraction.dat")

    gvec_chartmap = None
    cf_gvec = None
    if use_gvec:
        print(f"[{name}] Step 4: build GVEC chartmap from common VMEC input")
        gvec_chartmap = build_gvec_chartmap(
            wout,
            gvec_dir,
            minimize_tol=gvec_minimize_tol,
            max_iter=gvec_max_iter,
            total_iter=gvec_total_iter,
        )
        assert_boozer_chartmap_file(gvec_chartmap)
        shutil.copy2(start_boozer, gvec_dir / "start.dat")
        (gvec_dir / "simple.in").write_text(
            simple_in_chartmap("boozer_chartmap.nc", trace_time, facE_al),
            encoding="utf-8",
        )
        if list(gvec_dir.glob("wout*")):
            raise RuntimeError(f"{name}: GVEC run directory unexpectedly contains VMEC files")

        print(f"[{name}] Step 5: GVEC chartmap SIMPLE run")
        run_cmd([str(SIMPLE_X)], cwd=gvec_dir, label=f"{name} GVEC chartmap")
        cf_gvec = load_confined_fraction(gvec_dir / "confined_fraction.dat")

    print(f"[{name}] Step 6: field/orbit comparison artifacts")
    run_roundtrip(
        wout,
        chartmap_nc,
        "export",
        eq_dir / "roundtrip_vmec" / "vmec",
        f"{name}: VMEC direct vs VMEC chartmap",
        "VMEC chartmap",
        vmec_roundtrip_field_tol,
        "1.0e-6",
    )
    field_metrics: dict[str, float] = {}
    if use_gvec:
        run_roundtrip(
            wout,
            gvec_chartmap,
            "external",
            eq_dir / "roundtrip_gvec" / "gvec",
            f"{name}: VMEC direct vs GVEC chartmap",
            "GVEC chartmap",
            gvec_roundtrip_field_tol,
            "1.0e-6",
        )

        field_metrics = plot_field_component_comparison(
            eq_dir / f"{name.lower()}_field_components.png",
            name,
            chartmap_nc,
            gvec_chartmap,
            "VMEC chartmap",
            "GVEC chartmap",
        )
        plot_surface_comparison(
            eq_dir / f"{name.lower()}_surface_comparison.png",
            name,
            chartmap_nc,
            gvec_chartmap,
            "VMEC chartmap",
            "GVEC chartmap",
        )

    series = {
        "VMEC-Boozer": cf_vmec,
        "VMEC chartmap": cf_chartmap,
    }
    if use_gvec and cf_gvec is not None:
        series["GVEC chartmap"] = cf_gvec
    min_len = min(len(values) for values in series.values())
    metrics_by_label: dict[str, dict[str, np.ndarray | float]] = {}
    max_diffs = {"total": 0.0, "pass": 0.0, "trap": 0.0}
    ref_metrics = confined_metrics(series["VMEC-Boozer"][:min_len])
    tol = 1.0 / float(ref_metrics["npart"])
    for label, values in series.items():
        metrics_by_label[label] = confined_metrics(values[:min_len])

    labels = list(metrics_by_label)
    for idx, left in enumerate(labels):
        for right in labels[idx + 1:]:
            left_metrics = metrics_by_label[left]
            right_metrics = metrics_by_label[right]
            diff_total = float(
                np.max(np.abs(np.asarray(left_metrics["total_conf"]) - np.asarray(right_metrics["total_conf"])))
            )
            diff_pass = float(
                np.max(np.abs(np.asarray(left_metrics["pass_conf"]) - np.asarray(right_metrics["pass_conf"])))
            )
            diff_trap = float(
                np.max(np.abs(np.asarray(left_metrics["trap_conf"]) - np.asarray(right_metrics["trap_conf"])))
            )
            print(f"[{name}] max |total diff| {left} vs {right} = {diff_total:.6e}")
            print(f"[{name}] max |pass diff|  {left} vs {right} = {diff_pass:.6e}")
            print(f"[{name}] max |trap diff|  {left} vs {right} = {diff_trap:.6e}")
            max_diffs["total"] = max(max_diffs["total"], diff_total)
            max_diffs["pass"] = max(max_diffs["pass"], diff_pass)
            max_diffs["trap"] = max(max_diffs["trap"], diff_trap)

    plot_case_loss_comparison(eq_dir / f"e2e_{name}.png", name, metrics_by_label)

    summary = {
        label: summarize_confined_fraction(series[label][:min_len])
        for label in series
    }
    summary["executed_modes"] = list(series)
    summary["comparison"] = {
        "tol": tol,
        "max_total_diff": max_diffs["total"],
        "max_pass_diff": max_diffs["pass"],
        "max_trap_diff": max_diffs["trap"],
    }
    summary["field_metrics"] = field_metrics
    write_json(eq_dir / "loss_summary.json", summary)

    ok = all(value <= tol for value in max_diffs.values())
    if not ok:
        print(f"[{name}] FAIL: confined-fraction mismatch exceeds {tol:.6e}")
    else:
        print(f"[{name}] PASS")
    return ok, summary


def main() -> None:
    for path, label in [
        (SIMPLE_X, "simple.x"),
        (TOOL_X, "export tool"),
        (ROUNDTRIP_X, "roundtrip comparison tool"),
        (GVEC_TO_CHARTMAP, "GVEC conversion script"),
    ]:
        if not path.exists():
            raise SystemExit(f"Missing {label}: {path}")

    outdir = Path.cwd() / "boozer_chartmap_e2e"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    case_filter = os.environ.get("SIMPLE_BOOZER_CASES", "").strip()
    if case_filter:
        selected = {item.strip().upper() for item in case_filter.split(",") if item.strip()}
        equilibria = [eq for eq in EQUILIBRIA if str(eq["name"]).upper() in selected]
    else:
        equilibria = EQUILIBRIA

    all_ok = True
    combined_metrics: dict[str, dict[str, dict[str, np.ndarray | float]]] = {}
    summaries: dict[str, dict[str, dict[str, float]]] = {}
    for eq in equilibria:
        print(f"\n=== {eq['name']} ===")
        ok, summary = run_single_equilibrium(eq, outdir)
        series = {}
        eq_dir = outdir / str(eq["name"])
        for label, run_name in [
            ("VMEC-Boozer", "vmec_run"),
            ("VMEC chartmap", "vmec_chartmap_run"),
        ]:
            series[label] = confined_metrics(
                load_confined_fraction(eq_dir / run_name / "confined_fraction.dat")
            )
        if bool(eq.get("use_gvec", False)):
            series["GVEC chartmap"] = confined_metrics(
                load_confined_fraction(eq_dir / "gvec_chartmap_run" / "confined_fraction.dat")
            )
        combined_metrics[str(eq["name"])] = series
        summaries[str(eq["name"])] = summary
        all_ok = all_ok and ok

    plot_all_cases_total_losses(outdir / "all_e2e_losses.png", combined_metrics)
    write_json(outdir / "common_case_summaries.json", summaries)

    if not all_ok:
        raise SystemExit("One or more common-case end-to-end comparisons failed")
    print("ALL COMMON-CASE BOOZER CHARTMAP BENCHMARKS PASSED")


if __name__ == "__main__":
    main()
