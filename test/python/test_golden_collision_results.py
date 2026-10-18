"""Collision-only integrity contract and unchanged default/performance routing."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

SCRIPTS = Path(__file__).parents[1] / "golden_record"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location(
    "compare_collision_results", SCRIPTS / "compare_collision_results.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _case(path, current):
    path.mkdir(parents=True)
    starts = np.column_stack(
        [
            np.full(32, 0.3),
            np.arange(32) / 10,
            np.arange(32) / 20,
            np.ones(32),
            np.linspace(-0.8, 0.8, 32),
        ]
    )
    times = np.zeros((32, 10))
    times[:, 0] = np.arange(1, 33)
    times[:, 1] = 1e-4
    times[:, 2:5] = np.column_stack(
        [1 + np.arange(32) / 40, starts[:, 0], 1 - starts[:, 4] ** 2]
    )
    times[:, 5:10] = starts
    times[:, 8] = 0.99 + np.arange(32) * 1e-4 if current else 0.98
    exits = np.zeros((32, 5))
    exits[:, 0] = times[:, 0]
    exits[:, 2] = times[:, 1]
    # Same original marker loses at intentionally different stochastic times.
    times[0, 1] = 4e-5 if current else 5e-5
    times[0, 5] = 1.0
    exits[0, 1:3] = [1, times[0, 1]]
    np.savetxt(path / "start.dat", starts)
    np.savetxt(path / "times_lost.dat", times)
    np.savetxt(path / "orbit_exit_code.dat", exits)
    (path / "simple.in").write_text(
        "&config\nntestpart=32\ntrace_time=1d-4\nswcoll=.True.\ndeterministic=.True.\n/\n"
    )
    (path / "runtime_seconds.txt").write_text("2\n")


@pytest.fixture
def cases(tmp_path):
    ref, cur = tmp_path / "reference/collision", tmp_path / "current/collision"
    _case(ref, False)
    _case(cur, True)
    return ref, cur


def test_valid_decorrelated_collision_run(cases, capsys):
    assert MODULE.compare(*cases, 1.5e-7, 1e-12) == 0
    output = capsys.readouterr().out
    assert "not distribution validation" in output
    assert "Descriptive historical movement" in output


@pytest.mark.parametrize(
    "corruption",
    [
        "id",
        "initial",
        "start",
        "nan",
        "momentum",
        "pitch",
        "radius",
        "lcfs_accounting",
        "numerical",
        "skipped",
        "core",
        "wall",
        "fractional_code",
        "accounting",
        "premature",
        "negative_time",
        "late_time",
        "duplicate_kicks",
        "missing_row",
    ],
)
def test_current_corruption_fails(cases, corruption):
    ref, cur = cases
    times = np.loadtxt(cur / "times_lost.dat")
    exits = np.loadtxt(cur / "orbit_exit_code.dat")
    starts = np.loadtxt(cur / "start.dat")
    if corruption == "id":
        times[1, 0] = times[0, 0]
    elif corruption == "initial":
        times[1, 4] += 0.1
    elif corruption == "start":
        starts[1, 1] += 0.1
    elif corruption == "nan":
        times[1, 6] = np.nan
    elif corruption == "momentum":
        times[1, 8] = 0
    elif corruption == "pitch":
        times[1, 9] = 1.1
    elif corruption == "radius":
        times[1, 5] = 3
    elif corruption == "lcfs_accounting":
        times[0, 5] = 0.3
    elif corruption in ("numerical", "skipped", "core", "wall", "fractional_code"):
        exits[1, 1] = {
            "numerical": 102,
            "skipped": 3,
            "core": 4,
            "wall": 2,
            "fractional_code": 0.5,
        }[corruption]
    elif corruption == "accounting":
        exits[1, 2] /= 2
    elif corruption == "premature":
        times[1, 1] = exits[1, 2] = 5e-5
    elif corruption == "negative_time":
        times[0, 1] = exits[0, 2] = -1
    elif corruption == "late_time":
        times[0, 1] = exits[0, 2] = 2e-4
    elif corruption == "duplicate_kicks":
        times[1, 8] = times[2, 8]
    elif corruption == "missing_row":
        times = times[:-1]
    np.savetxt(cur / "start.dat", starts)
    np.savetxt(cur / "times_lost.dat", times)
    np.savetxt(cur / "orbit_exit_code.dat", exits)
    assert MODULE.compare(ref, cur, 1.5e-7, 1e-12) == 1


@pytest.mark.parametrize("which", [0, 1])
@pytest.mark.parametrize(
    "wall",
    [
        "'wall.stl'",
        "'/absolute/wall.stl'",
        "'wall!name.stl'",
        "'''quoted.stl'",
        "unquoted.stl",
    ],
)
def test_wall_configuration_fails_in_both_fixtures(cases, which, wall):
    path = cases[which] / "simple.in"
    path.write_text(path.read_text().replace("/\n", f"wall_input={wall}\n/\n"))
    assert MODULE.compare(*cases, 1.5e-7, 1e-12) == 1


def test_empty_wall_input_and_comment_do_not_enable_wall(cases):
    for case in cases:
        path = case / "simple.in"
        path.write_text(
            path.read_text().replace(
                "/\n", "wall_input='  '\nwall_units='m' ! wall_input='ignored.stl'\n/\n"
            )
        )
    assert MODULE.compare(*cases, 1.5e-7, 1e-12) == 0


@pytest.mark.parametrize("which", [0, 1])
def test_wall_character_substring_cannot_bypass_no_wall_gate(cases, which):
    path = cases[which] / "simple.in"
    path.write_text(path.read_text().replace("/\n", "wall_input(1:8)='wall.stl'\n/\n"))
    assert MODULE.compare(*cases, 1.5e-7, 1e-12) == 1


def _shell(ref, cur, enabled):
    env = dict(
        os.environ,
        GOLDEN_RECORD_STOCHASTIC_COLLISION=str(enabled),
        GOLDEN_RECORD_SKIP_CASES="",
        GOLDEN_RECORD_MAX_SLOWDOWN="1.2",
        GOLDEN_RECORD_MIN_REF_RUNTIME_S="1",
        GOLDEN_RECORD_MAX_ABS_SLOWDOWN_S=".2",
    )
    return subprocess.run(
        [
            "bash",
            str(SCRIPTS / "compare_golden_results.sh"),
            str(ref.parent),
            str(cur.parent),
            "collision",
        ],
        env=env,
        capture_output=True,
        text=True,
    )


def test_default_comparator_remains_strict_and_opt_in_is_explicit(cases):
    assert _shell(*cases, 0).returncode != 0
    assert _shell(*cases, 1).returncode == 0


def test_collision_performance_gate_remains_active(cases):
    ref, cur = cases
    (cur / "runtime_seconds.txt").write_text("3\n")
    result = _shell(ref, cur, 1)
    assert result.returncode != 0
    assert "performance regression" in result.stdout
