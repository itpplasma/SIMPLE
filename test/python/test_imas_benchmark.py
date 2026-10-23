"""Independent oracles for benchmark input validation and report decisions."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest


EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def load_example(name):
    spec = importlib.util.spec_from_file_location(name, EXAMPLES / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


benchmark = load_example("imas_benchmark")
regimes = load_example("imas_trajectory_regimes")


def test_fetch_checks_existing_cache(tmp_path, monkeypatch):
    path = tmp_path / "diii-d.geqdsk"
    path.write_bytes(b"corrupted equilibrium")
    monkeypatch.setattr(
        benchmark.urllib.request,
        "urlopen",
        lambda *a, **k: pytest.fail(
            "A cached input must be checked without downloading"
        ),
    )
    with pytest.raises(ValueError, match="Checksum mismatch"):
        benchmark.fetch("diii-d", tmp_path)


def test_revision_without_git_metadata(tmp_path):
    result = benchmark.revision(tmp_path)
    assert result == {"commit": None, "dirty": None}


@pytest.mark.parametrize(
    "option,value",
    [
        ("--energy", "nan"),
        ("--energy", "0"),
        ("--resolution", "7"),
        ("--perturbation", "0"),
        ("--horizon-tolerance", "inf"),
        ("--stat-markers", "-1"),
        ("--jobs", "0"),
    ],
)
def test_invalid_arguments_leave_no_output(tmp_path, option, value):
    output = tmp_path / "run"
    result = subprocess.run(
        [
            sys.executable,
            str(EXAMPLES / "imas_benchmark.py"),
            "--output",
            str(output),
            option,
            value,
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2, result.stderr
    assert not output.exists()


def orbit(t, pitch):
    t = np.asarray(t)
    return {
        "t": t,
        "r": np.full(len(t), 11.0),
        "z": np.zeros(len(t)),
        "phi": t.copy(),
        "pitch": np.asarray(pitch),
    }


def test_trapped_orbit_requires_a_complete_bounce():
    t = np.linspace(0, 1, 101)
    field = SimpleNamespace(axis=(10, 0), rho_pol=lambda r, z: r - 10)
    result = regimes.orbit_averages(orbit(t, np.cos(np.pi * t)), field)
    assert result["class"] == "incomplete"


def test_bounce_period_includes_exact_zero_samples():
    t = np.arange(17, dtype=float)
    pitch = np.tile([1.0, 0.0, -1.0, 0.0], 5)[:17]
    field = SimpleNamespace(axis=(10, 0), rho_pol=lambda r, z: r - 10)
    result = regimes.orbit_averages(orbit(t, pitch), field)
    assert result["class"] == "trapped"
    assert result["bounce_period_s"] == pytest.approx(4)


def test_incomplete_window_is_not_regular():
    record = {
        "unresolved": False,
        "lost": False,
        "classes": {"SIMPLE": "incomplete", "ASCOT5": "incomplete"},
        "t_pred_s": 1.0,
        "growth_exponent": None,
        "horizon_reached": {"SIMPLE": False, "ASCOT5": False},
    }
    assert regimes.marker_regime(record, 1.0) == "incomplete"


def test_horizon_crossing_at_last_sample_is_sensitive():
    record = {
        "unresolved": False,
        "lost": False,
        "classes": {"SIMPLE": "passing"},
        "min_abs_pitch": 0.5,
        "t_pred_s": 1.0,
        "growth_exponent": 3.0,
        "horizon_reached": {"SIMPLE": True, "ASCOT5": False},
    }
    assert regimes.marker_regime(record, 1.0) == "sensitive"


def test_unresolved_cohort_cannot_pass_statistics():
    end = {
        "resolved": np.zeros(4, dtype=bool),
        "lost": np.zeros(4, dtype=bool),
        **{key: np.zeros(4) for key in ("time", "r", "z", "pitch")},
    }
    field = SimpleNamespace(axis=(0, 0), rho_pol=lambda r, z: r)
    result = regimes.statistics(end, end, 1.0, field, ["SIMPLE", "ASCOT5"])
    assert result["markers_compared"] == 0
    assert result["accepted"] is None


def test_all_unresolved_orbits_still_produce_a_report(tmp_path, monkeypatch):
    empty = orbit([], [])
    empty.update(status="numerical(104)", loss_time=None, strike=None)
    nominal = orbit(np.linspace(0, 1, 9), np.ones(9))
    nominal.update(status="confined", loss_time=None, strike=None, speed=np.ones(9))
    calls = iter([nominal, nominal, empty, empty])
    monkeypatch.setattr(regimes, "simple_orbits", lambda *a, **k: ([next(calls)], 1.0))
    monkeypatch.setattr(
        regimes, "ascot_orbits", lambda eq, p, *a: ([empty] * len(p), {})
    )
    monkeypatch.setattr(
        regimes,
        "CommonField",
        lambda eq: SimpleNamespace(axis=(0, 0), rho_pol=lambda r, z: r),
    )
    result = regimes.compare(
        None,
        "unused",
        np.array([[1.0, 0.0, 0.0, 1e4, 0.5]]),
        2,
        1,
        duration=1,
        step=0.001,
        output=tmp_path,
    )
    assert result["summary"]["regime_counts"] == {"unresolved": 1}
    assert result["summary"]["t_pred_s_quantiles_resolved"] is None
    assert result["summary"]["max_start_position_gap_m"] is None
    regimes.write(result, tmp_path)


@pytest.mark.parametrize("version", ["3.42.2", "4.1.1"])
def test_analytic_equilibrium_cocos_and_units(version):
    pytest.importorskip("imas")
    from pysimple.imas.equilibrium import _geqdsk_input

    gfile = _geqdsk_input(benchmark.circular_equilibrium(version=version))
    r, z = np.meshgrid(np.linspace(8.7, 11.3, 129), np.linspace(-1.3, 1.3, 129))
    np.testing.assert_allclose(gfile["PsiVs"], 0.5 * ((r - 10) ** 2 + z**2), atol=1e-15)
    np.testing.assert_allclose(gfile["fprof"], 20)
    np.testing.assert_allclose(
        gfile["qprof"], 20 / np.sqrt(100 - np.linspace(0, 1, 129))
    )


def test_revision_does_not_attribute_an_enclosing_repository():
    assert benchmark.revision(EXAMPLES) == {"commit": None, "dirty": None}


def test_suite_retains_every_child_failure(tmp_path, monkeypatch):
    def child(command, **kwargs):
        return SimpleNamespace(returncode=7 if command[-1] == "1e-08" else 0)

    monkeypatch.setattr(benchmark.subprocess, "run", child)
    args = SimpleNamespace(suite="collisions", copies=1, jobs=2, cache=tmp_path)
    with pytest.raises(SystemExit) as error:
        benchmark.suite(args, tmp_path)
    assert error.value.code == 1
    report = json.loads((tmp_path / "suite.json").read_text())
    assert report["planned"] == len(report["runs"]) == 8
    assert sorted(r["exit_code"] for r in report["runs"]) == [0] * 4 + [7] * 4


@pytest.mark.parametrize("cocos,flux_scale", [(1, -1), (3, 1), (17, 2 * np.pi)])
@pytest.mark.parametrize("version", ["3.42.2", "4.1.1"])
def test_geqdsk_cocos_preserves_physical_field(monkeypatch, cocos, flux_scale, version):
    pytest.importorskip("imas")
    import libneo.eqdsk_base
    from pysimple.imas.equilibrium import _geqdsk_input

    expected = _geqdsk_input(benchmark.circular_equilibrium())
    source = dict(expected)
    source.update(R=np.linspace(8.7, 11.3, 129), Z=np.linspace(-1.3, 1.3, 129))
    for key in ("PsiVs", "PsiaxisVs", "PsiedgeVs"):
        source[key] = expected[key] * flux_scale
    monkeypatch.setattr(libneo.eqdsk_base, "read_eqdsk", lambda path: source)
    actual = _geqdsk_input(
        benchmark.equilibrium_from_geqdsk("unused", cocos, version=version)
    )
    for key in ("PsiVs", "PsiaxisVs", "PsiedgeVs", "fprof", "qprof", "Ip"):
        np.testing.assert_allclose(actual[key], expected[key], rtol=1e-12, atol=1e-14)


def test_deterministic_tolerance_uses_guiding_centre_path(tmp_path, monkeypatch):
    t = np.linspace(0, 1, 101)
    theta = 6 * np.pi * t
    native = orbit(t, np.full(len(t), 0.5))
    native.update(
        r=10 + 0.1 * np.cos(theta),
        z=0.1 * np.sin(theta),
        phi=np.zeros(len(t)),
        status="confined",
        loss_time=None,
        strike=None,
        speed=np.full(len(t), 1e6),
    )
    other = {**native, "r": native["r"] + 0.001}
    monkeypatch.setattr(regimes, "simple_orbits", lambda *a, **k: ([native], 1e6))
    monkeypatch.setattr(
        regimes, "ascot_orbits", lambda eq, p, *a: ([other] * len(p), {})
    )
    monkeypatch.setattr(
        regimes,
        "CommonField",
        lambda eq: SimpleNamespace(
            axis=(10, 0), rho_pol=lambda r, z: np.hypot(r - 10, z)
        ),
    )
    monkeypatch.setattr(regimes, "invariant_drift", lambda *a: {})
    report = regimes.compare(
        None,
        "unused",
        np.array([[10.1, 0.0, 0.0, 1e4, 0.5]]),
        2,
        1,
        duration=1,
        step=0.001,
        output=tmp_path,
    )
    marker = report["markers"][0]
    assert marker["regime"] == "regular_passing"
    assert not marker["deterministic_accepted"]
    window = marker["pointwise"]["relativistic_equivalent"]
    assert window["path_length_m"] == pytest.approx(0.6 * np.pi, rel=0.002)
