#!/usr/bin/env python3
"""Driver coverage for standalone fast classification.

fast_class documents itself as "quit immediately after fast classification".
Standalone (tcut <= 0, class_plot = .False.) it is the only configuration that
classifies orbits without the Minkowski fractal cut, which fires at
kt == ntcut, and without disabling the early exit, which requires
.not. class_plot.

The oracle is the classification output itself: a classifying run owes a
class_parts.dat with one row per test particle, and a non-classifying run owes
none. Both expectations are independent of how the dispatch predicate is
spelled. A restart that keeps every finished marker owes the same
class_parts.dat and class_scores.dat as the run it resumes.

The early exit itself is checked against the same markers traced with the
Minkowski cut at trace_time instead: an orbit that fast classification stops
as J_parallel-regular must end elsewhere than the full trace, and every other
orbit must end where the full trace ends.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from simple_driver import assert_file_absent, assert_row_count, run_simple_case

NTESTPART = 4

COMMON = """
&config
multharm = 3
contr_pp = -1d10
trace_time = 1d-4
ntestpart = 4
nper = 1
npoiper = 4
npoiper2 = 32
ntimstep = 4
num_surf = 1
sbeg = 0.3
notrace_passing = 1
netcdffile = 'wout.nc'
isw_field_type = 2
integmode = 3
deterministic = .True.
/
"""

# (name, extra namelist lines, classification expected)
CASES = (
    ("fast_class_standalone", "fast_class = .True.\ntcut = -1d0", True),
    ("fast_class_with_cut", "fast_class = .True.\ntcut = 5d-5", True),
    ("fast_class_with_plot", "fast_class = .True.\nclass_plot = .True.", True),
    ("fast_class_absent", "tcut = -1d0", False),
)


def config_with(additions: str) -> str:
    return COMMON.replace(
        "deterministic = .True.",
        f"deterministic = .True.\n{additions}",
    )


# Marker 4 at s = 0.6 is trapped and J_parallel-regular within this trace.
EARLY_STOP = """
&config
multharm = 3
notrace_passing = 1
ntestpart = {ntestpart}
sbeg = 0.6
npoiper2 = 128
netcdffile = 'wout.nc'
isw_field_type = 2
trace_time = 1.5d-2
deterministic = .True.
{extra}
/
"""


def check_early_stop(simple_x: Path, wout: Path) -> None:
    fast = run_simple_case(
        simple_x, wout, "fast_class_early_stop",
        EARLY_STOP.format(ntestpart=NTESTPART,
                          extra="fast_class = .True.\ntcut = -1d0"), timeout=300)
    full = run_simple_case(
        simple_x, wout, "fast_class_full_trace",
        EARLY_STOP.format(ntestpart=NTESTPART, extra="tcut = 1.5d-2"), timeout=300)

    jpar = [row[3] for row in table(fast / "class_parts.dat")]
    exit_code = [row[1] for row in table(fast / "orbit_exit_code.dat")]
    fast_end = [row[5:] for row in table(fast / "times_lost.dat")]
    full_end = [row[5:] for row in table(full / "times_lost.dat")]

    stopped = [i + 1 for i in range(NTESTPART)
               if jpar[i] == "1" and exit_code[i] == "0"]
    if not stopped:
        raise AssertionError("no J_parallel-regular confined marker to stop early")
    moved = [i + 1 for i in range(NTESTPART) if fast_end[i] != full_end[i]]
    if moved != stopped:
        raise AssertionError(
            f"orbits ending off the full trace {moved}, "
            f"expected the early-stopped {stopped}")


def table(path: Path) -> list[list[str]]:
    return [line.split() for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def check_restart_keeps_classes(simple_x: Path, wout: Path) -> None:
    case_dir = run_simple_case(simple_x, wout, "fast_class_restart",
                               config_with("fast_class = .True.\ntcut = 5d-5"))
    first = {name: (case_dir / name).read_text(encoding="utf-8")
             for name in ("class_parts.dat", "class_scores.dat")}
    if all(line.split()[-1] == "0" for line in first["class_parts.dat"].splitlines()):
        raise AssertionError("fast_class_restart: no classified marker to restore")
    (case_dir / "simple.in").write_text(
        config_with("fast_class = .True.\ntcut = 5d-5\nrestart = .True."),
        encoding="utf-8")
    result = subprocess.run([str(simple_x)], cwd=case_dir, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=90, check=False)
    if result.returncode != 0:
        raise AssertionError(f"fast_class_restart: restart failed\n{result.stdout}")
    for name, text in first.items():
        if (case_dir / name).read_text(encoding="utf-8") != text:
            raise AssertionError(f"fast_class_restart: restart changed {name}")


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_fast_class_driver.py SIMPLE_X WOUT_FILE", file=sys.stderr)
        return 2

    simple_x = Path(sys.argv[1]).resolve()
    wout = Path(sys.argv[2]).resolve()

    for name, additions, classifies in CASES:
        case_dir = run_simple_case(simple_x, wout, name, config_with(additions))
        if classifies:
            assert_row_count(case_dir, "class_parts.dat", NTESTPART)
        else:
            assert_file_absent(case_dir, "class_parts.dat")

    check_early_stop(simple_x, wout)
    check_restart_keeps_classes(simple_x, wout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
