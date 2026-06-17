"""Tests for tools/step_convergence.py (statistical step-size convergence check)."""
import importlib.util
from pathlib import Path

import numpy as np

_TOOL = Path(__file__).resolve().parents[2] / "tools" / "step_convergence.py"
_spec = importlib.util.spec_from_file_location("step_convergence", _TOOL)
sc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sc)


def _write(path, loss_times):
    # times_lost.dat columns: index, loss_time, then padding columns.
    n = len(loss_times)
    arr = np.zeros((n, 10))
    arr[:, 0] = np.arange(1, n + 1)
    arr[:, 1] = loss_times
    np.savetxt(path, arr)


def test_identical_runs_converged(tmp_path):
    rng = np.random.default_rng(0)
    # half confined (-1), half lost at random positive times
    t = np.where(rng.random(1000) < 0.5, -1.0, rng.uniform(0.01, 1.0, 1000))
    _write(tmp_path / "c.dat", t)
    _write(tmp_path / "f.dat", t)
    r = sc.compare(tmp_path / "c.dat", tmp_path / "f.dat")
    assert r["bias"] == 0.0
    assert r["flip_rate"] == 0.0
    assert r["ks"] == 0.0
    assert r["converged"]


def test_known_flips_counted(tmp_path):
    n = 1000
    tc = np.full(n, -1.0)            # all confined at coarse
    tf = np.full(n, -1.0)
    # 7 particles confined@coarse become lost@fine; 3 the other way.
    tf[:7] = 0.5                      # conf -> lost  (n_cf = 7)
    tc[7:10] = 0.5                    # lost -> conf  (n_fc = 3)
    _write(tmp_path / "c.dat", tc)
    _write(tmp_path / "f.dat", tf)
    r = sc.compare(tmp_path / "c.dat", tmp_path / "f.dat")
    assert r["n_cf"] == 7
    assert r["n_fc"] == 3
    assert r["bias"] == (7 - 3) / n
    assert r["flip_rate"] == 10 / n


def test_large_bias_not_converged(tmp_path):
    n = 1000
    tc = np.full(n, -1.0)            # all confined at coarse
    tf = np.where(np.arange(n) < 200, 0.3, -1.0)  # 200 lost at fine
    _write(tmp_path / "c.dat", tc)
    _write(tmp_path / "f.dat", tf)
    r = sc.compare(tmp_path / "c.dat", tmp_path / "f.dat")
    assert abs(r["bias"]) > r["mc_err"]
    assert not r["converged"]
