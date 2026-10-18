"""Independent oracles for report statistics, fits, and failure handling."""

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.stats import binomtest

_SPEC = importlib.util.spec_from_file_location(
    "verification_report",
    Path(__file__).resolve().parents[2] / "tools/generate_verification_report.py",
)
vr = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(vr)


@pytest.mark.parametrize(
    "k,n", [(0, 1), (1, 1), (0, 100), (100, 100), (7, 32), (128, 256)]
)
def test_wilson_matches_scipy(k, n):
    reference = binomtest(k, n).proportion_ci(method="wilson")
    np.testing.assert_allclose(
        vr.wilson_interval(k, n), [reference.low, reference.high], atol=2e-15
    )


@pytest.mark.parametrize("k,n", [(-1, 10), (11, 10), (0, 0)])
def test_invalid_counts(k, n):
    with pytest.raises(ValueError):
        vr.wilson_interval(k, n)


def test_fit_known_order():
    h = np.array([0.1, 0.05, 0.025, 0.0125])
    order, _ = vr.fit_order(h, 3 * h**2)
    assert order == pytest.approx(2.0)


@pytest.mark.parametrize("errors", [[0, 1, 2], [1, np.nan, 2], [1, np.inf, 2]])
def test_nonfinite_or_zero_error_cannot_pass(errors):
    with pytest.raises(ValueError):
        vr.fit_order([1, 0.5, 0.25], errors)


def test_resolved_loss_outcomes():
    np.testing.assert_array_equal(
        vr.loss_mask([0, 1, 2], [1, 0.3, 1], 1), [False, True, True]
    )


@pytest.mark.parametrize(
    "codes,times",
    [
        ([101], [np.nan]),
        ([3], [-1]),
        ([0], [0.5]),
        ([1], [2]),
        ([1], [0]),
        ([1.5], [0.3]),
        ([], []),
    ],
)
def test_unresolved_or_inconsistent_is_not_binomial(codes, times):
    with pytest.raises(ValueError):
        vr.loss_mask(codes, times, 1)


def test_missing_metric_is_failure():
    with pytest.raises(ValueError, match="missing measured metric"):
        vr.metric("some unrelated test passed", "required error:")


def test_timeout_preserves_log_and_fails(tmp_path):
    runner = vr.Runner(tmp_path, 0.5)
    with pytest.raises(RuntimeError, match="time budget"):
        runner.run(
            [
                sys.executable,
                "-c",
                "import time; print('started', flush=True); time.sleep(10)",
            ],
            tmp_path,
            tmp_path / "timeout.log",
        )
    assert "started" in (tmp_path / "timeout.log").read_text()


def test_empty_or_partial_test_selection_fails(tmp_path):
    class MissingTests:
        output = tmp_path

        def run(self, *args):
            return "No tests were found!!!"

    with pytest.raises(ValueError, match="did not pass or was not registered"):
        vr.consumer_checks(MissingTests())


@pytest.mark.parametrize("exception", [False, True])
def test_report_failure_returns_nonzero_and_preserves_explanation(
    monkeypatch, tmp_path, exception
):
    output = tmp_path / "build/report"
    monkeypatch.setattr(vr, "REPO", tmp_path)
    monkeypatch.setattr(sys, "argv", ["report", "--output-dir", str(output)])
    monkeypatch.setattr(vr, "provenance", lambda: [])

    def checks(runner):
        if exception:
            raise ValueError("missing numerical evidence")
        return [
            vr.row("deliberate numerical failure", "oracle", "bound", "error", False)
        ]

    monkeypatch.setattr(vr, "consumer_checks", checks)
    monkeypatch.setattr(vr, "convergence", lambda *args: [])
    monkeypatch.setattr(vr, "uncertainty", lambda *args: [])
    assert vr.main() == 1
    report = (output / "verification-report.md").read_text()
    assert "FAIL" in report
    assert (
        "missing numerical evidence" if exception else "deliberate numerical failure"
    ) in report
