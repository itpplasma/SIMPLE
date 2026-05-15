#!/usr/bin/env python3
"""Targeted driver coverage for bminmax cache lifecycle."""

from __future__ import annotations

import sys
from pathlib import Path

from simple_driver import (
    assert_file_absent,
    assert_log_order,
    assert_row_count,
    log_text,
    run_simple_case,
)

FIELD_INIT_MARKER = "Field type initialization completed"
BMINMAX_INIT_MARKER = "Bmin/Bmax initialization completed"

CLASSIFIER_CONFIG = """\
class_plot = .True.
fast_class = .True.
notrace_passing = 1
tcut = -1d0
cut_in_per = 0d0
"""


def assert_no_bminmax_cache(case_dir: Path) -> None:
    assert_file_absent(case_dir, "bminmax.dat")
    if BMINMAX_INIT_MARKER in log_text(case_dir):
        raise AssertionError(
            f"{case_dir.name}: bminmax cache was initialized unexpectedly"
        )


def config_with(base: str, additions: str) -> str:
    if not additions:
        return base
    return base.replace(
        "deterministic = .True.",
        f"deterministic = .True.\n{additions}",
    )


def assert_expected_cache(case_dir: Path, expect_cache: bool) -> None:
    if expect_cache:
        assert_log_order(case_dir, FIELD_INIT_MARKER, BMINMAX_INIT_MARKER)
    else:
        assert_no_bminmax_cache(case_dir)


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

    cases = [
        ("bminmax_single_surface", "", False, 0),
        ("bminmax_volume", "startmode = 5\nnum_surf = 0", True, 0),
        (
            "bminmax_classifier_num_surf_zero",
            f"num_surf = 0\n{CLASSIFIER_CONFIG}",
            False,
            0,
        ),
        (
            "bminmax_classifier_multi_surface",
            f"num_surf = 2\n{CLASSIFIER_CONFIG}",
            True,
            101,
        ),
    ]

    for name, additions, expect_cache, expected_rows in cases:
        case_dir = run_simple_case(
            simple_x,
            wout,
            name,
            config_with(common, additions),
        )
        assert_expected_cache(case_dir, expect_cache)
        if expected_rows:
            assert_row_count(case_dir, "bminmax.dat", expected_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
