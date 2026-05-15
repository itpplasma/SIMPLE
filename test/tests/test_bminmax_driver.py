#!/usr/bin/env python3
"""Targeted driver coverage for bminmax cache lifecycle."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def run_case(simple_x: Path, wout: Path, name: str, config: str) -> Path:
    work_root = Path(os.environ.get("CTEST_BINARY_DIRECTORY", Path.cwd()))
    case_dir = Path(tempfile.mkdtemp(prefix=f"{name}_", dir=work_root))
    (case_dir / "simple.in").write_text(config, encoding="utf-8")
    try:
        (case_dir / "wout.nc").symlink_to(wout)
    except OSError:
        shutil.copy2(wout, case_dir / "wout.nc")

    result = subprocess.run(
        [str(simple_x)],
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=90,
        check=False,
    )
    (case_dir / "simple.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        raise AssertionError(
            f"{name}: simple.x failed with {result.returncode}\n{result.stdout}"
        )
    return case_dir


def assert_log_order(case_dir: Path, before: str, after: str) -> None:
    log = (case_dir / "simple.log").read_text(encoding="utf-8")
    before_idx = log.find(before)
    after_idx = log.find(after)
    if before_idx < 0 or after_idx < 0:
        raise AssertionError(f"{case_dir.name}: missing log markers\n{log}")
    if before_idx > after_idx:
        raise AssertionError(
            f"{case_dir.name}: wrong log order for {before!r} and {after!r}"
        )


def assert_has_bminmax(case_dir: Path) -> None:
    bminmax = case_dir / "bminmax.dat"
    if not bminmax.exists():
        raise AssertionError(f"{case_dir.name}: expected bminmax.dat")
    rows = [
        line
        for line in bminmax.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(rows) != 101:
        raise AssertionError(
            f"{case_dir.name}: expected 101 bminmax rows, got {len(rows)}"
        )


def assert_no_bminmax(case_dir: Path) -> None:
    if (case_dir / "bminmax.dat").exists():
        raise AssertionError(f"{case_dir.name}: bminmax.dat should not be written")
    log = (case_dir / "simple.log").read_text(encoding="utf-8")
    if "Bmin/Bmax initialization completed" in log:
        raise AssertionError(
            f"{case_dir.name}: bminmax cache was initialized unexpectedly"
        )


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_bminmax_driver.py SIMPLE_X WOUT_FILE", file=sys.stderr)
        return 2

    simple_x = Path(sys.argv[1]).resolve()
    wout = Path(sys.argv[2]).resolve()

    common = """
&config
multharm = 3
contr_pp = -1d10
trace_time = 1d-5
ntestpart = 2
nper = 1
npoiper = 4
npoiper2 = 32
ntimstep = 2
sbeg = 0.2, 0.4
netcdffile = 'wout.nc'
isw_field_type = 2
deterministic = .True.
/
"""

    single_surface = run_case(simple_x, wout, "bminmax_single_surface", common)
    assert_no_bminmax(single_surface)

    volume = run_case(
        simple_x,
        wout,
        "bminmax_volume",
        common.replace(
            "deterministic = .True.",
            "deterministic = .True.\nstartmode = 5\nnum_surf = 0",
        ),
    )
    assert_log_order(
        volume,
        "Field type initialization completed",
        "Bmin/Bmax initialization completed",
    )

    classifier_num_surf_zero = run_case(
        simple_x,
        wout,
        "bminmax_classifier_num_surf_zero",
        common.replace(
            "deterministic = .True.",
            "deterministic = .True.\nnum_surf = 0\nclass_plot = .True.\n"
            "fast_class = .True.\n"
            "notrace_passing = 1\ntcut = -1d0\ncut_in_per = 0d0",
        ),
    )
    assert_no_bminmax(classifier_num_surf_zero)

    classifier_multi_surface = run_case(
        simple_x,
        wout,
        "bminmax_classifier_multi_surface",
        common.replace(
            "deterministic = .True.",
            "deterministic = .True.\nnum_surf = 2\nclass_plot = .True.\n"
            "fast_class = .True.\n"
            "notrace_passing = 1\ntcut = -1d0\ncut_in_per = 0d0",
        ),
    )
    assert_log_order(
        classifier_multi_surface,
        "Field type initialization completed",
        "Bmin/Bmax initialization completed",
    )
    assert_has_bminmax(classifier_multi_surface)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
