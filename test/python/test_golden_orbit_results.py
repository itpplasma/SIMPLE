"""Tests for accounting-aware golden orbit comparison."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).parents[1] / "golden_record" / "compare_orbit_results.py"
SPEC = importlib.util.spec_from_file_location("compare_orbit_results", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _write_case(path: Path, times: np.ndarray, exits: np.ndarray) -> None:
    path.mkdir()
    np.savetxt(path / "times_lost.dat", times)
    np.savetxt(path / "orbit_exit_code.dat", exits)
    (path / "simple.in").write_text(
        "&config\n  trace_time = 1d-4\n/\n", encoding="utf-8"
    )


def _base_tables() -> tuple[np.ndarray, np.ndarray]:
    times = np.array(
        [
            [1.0, 1.0e-4, 1.0, 0.3, 0.2, 0.4, 0.1, 0.2, 1.0, 0.0],
            [2.0, 4.0e-5, -1.0, 0.3, 0.2, 1.0, 0.2, 0.3, 1.0, 0.0],
        ]
    )
    exits = np.array(
        [
            [1.0, 0.0, 1.0e-4, 0.0, 0.0],
            [2.0, 1.0, 4.0e-5, 0.0, 0.0],
        ]
    )
    return times, exits


def test_exact_orbit_results_match(tmp_path: Path) -> None:
    times, exits = _base_tables()
    _write_case(tmp_path / "ref", times, exits)
    _write_case(tmp_path / "cur", times, exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 0


def test_numerical_exit_may_recover_to_survivor(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    ref_times[0, 1] = np.nan
    ref_exits[0, 1:3] = [104.0, np.nan]
    cur_times, cur_exits = _base_tables()
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 3


def test_numerical_exit_may_recover_to_physical_loss(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    ref_times[0, 1] = np.nan
    ref_exits[0, 1:3] = [103.0, np.nan]
    cur_times, cur_exits = _base_tables()
    cur_times[0, 1] = 7.0e-5
    cur_exits[0, 1:3] = [2.0, 7.0e-5]
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 3


def test_premature_recovered_completion_fails(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    ref_times[0, 1] = np.nan
    ref_exits[0, 1:3] = [104.0, np.nan]
    cur_times, cur_exits = _base_tables()
    cur_times[0, 1] = 8.0e-5
    cur_exits[0, 2] = 8.0e-5
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 1


def test_ordinary_physics_drift_still_fails(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    cur_times, cur_exits = _base_tables()
    cur_times[1, 1] = 5.0e-5
    cur_exits[1, 2] = 5.0e-5
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 1


def test_valid_survivor_endpoint_may_decorrelate(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    cur_times, cur_exits = _base_tables()
    cur_times[0, 5:10] = [0.8, 0.2, 0.4, 0.9, -0.1]
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 0


def test_invalid_survivor_endpoint_fails(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    cur_times, cur_exits = _base_tables()
    cur_times[0, 5] = np.nan
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 1


def test_unrelated_columns_on_recovered_marker_still_fail(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    ref_times[0, 1] = np.nan
    ref_exits[0, 1:3] = [105.0, np.nan]
    cur_times, cur_exits = _base_tables()
    cur_times[0, 3] = 0.4
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 1


def test_recovered_marker_may_have_a_new_final_state(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    ref_times[0, 1] = np.nan
    ref_times[0, 5:10] = np.nan
    ref_exits[0, 1:3] = [105.0, np.nan]
    cur_times, cur_exits = _base_tables()
    cur_times[0, 5:10] = [0.8, 0.2, 0.4, 0.9, -0.1]
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 3


def test_non_numerical_reference_transition_fails(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    ref_times[0, 1] = 6.0e-5
    ref_exits[0, 1:3] = [2.0, 6.0e-5]
    cur_times, cur_exits = _base_tables()
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 1


def test_invalid_physical_endpoint_may_be_corrected(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    ref_times[1, 5] = 5.0e5
    cur_times, cur_exits = _base_tables()
    cur_times[1, 1] = 1.0e-4
    cur_times[1, 5:10] = [0.8, 0.2, 0.4, 0.9, -0.1]
    cur_exits[1, 1:3] = [0.0, 1.0e-4]
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 3


def test_invalid_physical_endpoint_requires_valid_correction(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    ref_times[1, 5] = 5.0e5
    cur_times, cur_exits = _base_tables()
    cur_times[1, 1] = 1.0e-4
    cur_times[1, 5] = np.nan
    cur_exits[1, 1:3] = [0.0, 1.0e-4]
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 1


def test_near_lcfs_survivor_may_become_lost(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    ref_times[0, 5] = 0.9999
    cur_times, cur_exits = _base_tables()
    cur_times[0, 1] = 2.0e-5
    cur_times[0, 5] = 1.0
    cur_exits[0, 1:3] = [1.0, 2.0e-5]
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 3


def test_interior_survivor_may_not_become_lost(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    cur_times, cur_exits = _base_tables()
    cur_times[0, 1] = 2.0e-5
    cur_times[0, 5] = 1.0
    cur_exits[0, 1:3] = [1.0, 2.0e-5]
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 1


def test_duplicate_loss_time_disagreement_fails(tmp_path: Path) -> None:
    ref_times, ref_exits = _base_tables()
    cur_times, cur_exits = _base_tables()
    cur_exits[0, 2] = 9.0e-5
    _write_case(tmp_path / "ref", ref_times, ref_exits)
    _write_case(tmp_path / "cur", cur_times, cur_exits)
    assert MODULE.compare(tmp_path / "ref", tmp_path / "cur", 1.0e-7, 1.0e-12) == 1
