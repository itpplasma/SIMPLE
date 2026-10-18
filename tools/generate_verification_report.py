#!/usr/bin/env python3
"""Run bounded checks of built SIMPLE; write Markdown, raw data, logs and PNGs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np

REPO = Path(__file__).resolve().parents[1]
TESTS = (
    "generate_test_boozer_chartmap_data",
    "test_profiles",
    "test_orbit_invariants",
    "test_boozer_chartmap",
    "test_chartmap_scaling",
    "test_boozer_chartmap_roundtrip",
)
NIST = "https://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm"


def wilson_interval(k: int, n: int) -> tuple[float, float]:
    """Two-sided 95% Wilson score interval, including k=0 and k=n."""
    if n <= 0 or not 0 <= k <= n:
        raise ValueError("require n > 0 and 0 <= k <= n")
    z = 1.959963984540054
    p = k / n
    denominator = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denominator
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return max(0.0, centre - half), min(1.0, centre + half)


def fit_order(steps, errors) -> tuple[float, float]:
    steps, errors = np.asarray(steps), np.asarray(errors)
    if (
        steps.ndim != 1
        or steps.shape != errors.shape
        or len(steps) < 3
        or not np.isfinite([steps, errors]).all()
        or np.any(steps <= 0)
        or np.any(errors <= 0)
        or len(np.unique(steps)) != len(steps)
    ):
        raise ValueError(
            "order fit needs at least three distinct positive finite steps/errors"
        )
    slope, intercept = np.polyfit(np.log(steps), np.log(errors), 1)
    return float(slope), float(intercept)


def loss_mask(codes, times, trace_time: float) -> np.ndarray:
    """Accept only resolved, traced markers; never turn numerical failures into losses."""
    codes, times = np.asarray(codes), np.asarray(times)
    if (
        codes.shape != times.shape
        or codes.size == 0
        or not np.isfinite(times).all()
        or not np.isin(codes, [0, 1, 2]).all()
        or np.any(times <= 0)
        or np.any(times > trace_time * (1 + 1e-9))
        or not np.allclose(times[codes == 0], trace_time, rtol=1e-9, atol=0)
    ):
        raise ValueError(
            "unresolved, skipped or inconsistent marker outcomes; no binomial estimate"
        )
    return np.isin(codes, [1, 2])


class Runner:
    def __init__(self, output: Path, seconds: float):
        self.output = output
        self.deadline = time.monotonic() + seconds

    def run(self, command, cwd: Path, log: Path) -> str:
        remaining = self.deadline - time.monotonic()
        if remaining <= 0:
            raise RuntimeError("report time budget exhausted")
        env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
        # Kill the whole make/CTest process group on timeout, preserving partial logs.
        with log.open("w") as stream:
            stream.write(f"$ {shlex.join(map(str, command))}\n")
            stream.flush()
            with subprocess.Popen(
                command,
                cwd=cwd,
                env=env,
                stdout=stream,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            ) as proc:
                try:
                    code = proc.wait(timeout=remaining)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                    raise RuntimeError(f"time budget exceeded; see {log}") from None
        if code:
            raise RuntimeError(f"command exited {code}; see {log}")
        return log.read_text()


def row(name, oracle, tolerance, measured, passed=True):
    return dict(
        check=name,
        oracle=oracle,
        tolerance=tolerance,
        measured=measured,
        status="PASS" if passed else "FAIL",
    )


def metric(log: str, pattern: str) -> float:
    match = re.search(pattern + r"\s*([+-]?[\d.]+(?:[EeDd][+-]?\d+)?)", log)
    if match is None:
        raise ValueError(f"missing measured metric in consumer-tests.log: {pattern}")
    value = float(match[1].replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError(f"nonfinite consumer metric: {pattern}")
    return value


def consumer_checks(runner: Runner) -> list[dict]:
    # Make expands command-line variables; $$ preserves the CTest end anchor.
    pattern = "^(" + "|".join(TESTS) + ")$$"
    log = runner.run(
        ["make", "test", "CONFIG=Fast", f"TEST={pattern}", "VERBOSE=1"],
        REPO,
        runner.output / "consumer-tests.log",
    )
    for name in TESTS:
        if not re.search(r"Test\s+#\d+: " + re.escape(name) + r"\s+\.+\s+Passed", log):
            raise ValueError(
                f"{name} did not pass or was not registered; see consumer-tests.log"
            )
    ratio = metric(log, r"<p\^2> final/eq =")
    shape = metric(log, r"distribution stays Maxwellian, max reldev")
    slowing = metric(log, r"distribution matches p\^2/\(p\^3\+pc\^3\), max reldev")
    field = metric(log, r"max relative error Bmod:")
    orbit = metric(log, r"max \|s_direct - s_chartmap\|:")
    return [
        row(
            "Maxwellian collision fixed point",
            "<p²>=3/(2 enrat); Maxwellian shape; 8000 markers, 4000 steps",
            "0.95 < moment ratio < 1.05; populated-bin shape deviation < 0.1",
            f"ratio={ratio:.8g}; max shape deviation={shape:.6g}",
            0.95 < ratio < 1.05 and shape < 0.1,
        ),
        row(
            "Slowing-down distribution",
            "normalized p²/(p³+pc³) residence histogram",
            "max relative bin error <= 0.04",
            f"max relative error={slowing:.6g}",
            slowing <= 0.04,
        ),
        row(
            "Invariant reconstruction",
            "forward/inverse H0, J_perp, p_phi in analytic TEST field",
            "scaled residuals <= 1e-11, 1e-10, 1e-9 respectively",
            "all assertions passed; individual residuals not printed (bounds only)",
        ),
        row(
            "Chartmap coordinates and scaling",
            "rho→s→rho; B×5, length×3, flux×45",
            "roundtrip <= 1e-14; scaling 1e-12, B interpolation 5e-3, Cartesian inverse 1e-10",
            "all assertions passed; individual residuals not printed (bounds only)",
        ),
        row(
            "VMEC/Boozer export/reimport",
            "direct Boozer field and 500-step midpoint orbit; 50 field points",
            "max relative B error <= 8e-5; max absolute s difference <= 2e-9",
            f"B error={field:.6g}; s difference={orbit:.6g}",
            field <= 8e-5 and orbit <= 2e-9,
        ),
    ]


def simple_run(runner: Runner, name: str, settings: dict, starts=None):
    from netCDF4 import Dataset

    work = runner.output / name
    work.mkdir()
    config = dict(
        deterministic=".true.",
        contr_pp="-1d10",
        relerr="1d-12",
        output_orbits_macrostep=".true.",
        **settings,
    )
    (work / "simple.in").write_text(
        "&config\n"
        + "\n".join(f"{key} = {value}" for key, value in config.items())
        + "\n/\n"
    )
    if starts is not None:
        np.savetxt(work / "start.dat", starts, fmt="%.17e")
    runner.run([str(REPO / "build/simple.x"), "simple.in"], work, work / "run.log")
    outcomes = np.loadtxt(work / "orbit_exit_code.dat", ndmin=2)
    loss_mask(outcomes[:, 1], outcomes[:, 2], float(config["trace_time"]))
    with Dataset(work / "orbits.nc") as ds:
        trajectory = np.stack(
            [
                np.ma.filled(ds[key][:], np.nan)
                for key in ("s", "theta", "phi", "p_abs", "v_par")
            ]
        )
        times = np.ma.filled(ds["time"][:], np.nan)[:, 0]
    return trajectory, times, outcomes


def convergence(runner: Runner, plt) -> list[dict]:
    initial = None
    states, errors, steps = [], [], []
    traces = []
    for mode, count in [(5, 128), (5, 256), (3, 8), (3, 16), (3, 32), (3, 64)]:
        name = f"orbit-m{mode}-n{count}"
        trajectory, times, outcomes = simple_run(
            runner,
            name,
            dict(
                isw_field_type=-1,
                integmode=mode,
                trace_time=0.1,
                ntimstep=count + 1,
                npoiper2=1,
                ntestpart=1,
                sbeg=0.4,
                ran_seed=42,
            ),
        )
        start = np.loadtxt(runner.output / name / "start.dat")
        if initial is None:
            initial = start
        if (
            not np.array_equal(start, initial)
            or not np.isfinite(trajectory).all()
            or trajectory.shape != (5, count + 1, 1)
            or not np.allclose(
                times, np.linspace(0, 0.1, count + 1), atol=1e-14, rtol=1e-12
            )
            or outcomes[0, 1] != 0
        ):
            raise ValueError(
                f"incomplete or inconsistent convergence trajectory: {name}"
            )
        states.append(trajectory[:, -1, 0])
        traces.append(trajectory[:, :, 0])
    reference = states[1]
    # Gauss returns the last-stage radius (orbit_symplectic.f90), while
    # theta/phi are endpoint updates. Compare endpoint angles, never mix times.
    reference_gap = float(np.max(np.abs(states[0][1:3] - reference[1:3])))
    for count, state in zip((8, 16, 32, 64), states[2:]):
        steps.append(0.1 / count)
        errors.append(float(np.max(np.abs(state[1:3] - reference[1:3]))))
    slope, intercept = fit_order(steps, errors)
    energy = max(
        float(np.max(np.abs(t[3] ** 2 - initial[3] ** 2)) / initial[3] ** 2)
        for t in traces[2:]
    )
    magnetic_moment_error, momentum_error = 0.0, 0.0
    for r, theta, _, speed, pitch in traces[2:]:
        # Independent analytic definitions from field_can_test, B0=R0=ro0=1,
        # a=0.5. Axisymmetry conserves native toroidal canonical momentum.
        j_perp = speed**2 * (1 - pitch**2) / (1 - r * np.cos(theta))
        p_phi = -(r**2) / 2 + r**4 + speed * pitch * (1 + r * np.cos(theta))
        magnetic_moment_error = max(
            magnetic_moment_error,
            float(np.max(np.abs(j_perp - j_perp[0])) / abs(j_perp[0])),
        )
        momentum_error = max(momentum_error, float(np.max(np.abs(p_phi - p_phi[0]))))
    np.savetxt(
        runner.output / "convergence.csv",
        np.column_stack([steps, errors]),
        delimiter=",",
        header="normalized_step,max_absolute_angle_error_rad",
        comments="",
    )
    fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
    ax.loglog(steps, errors, "o", color="#0072B2", label="SIMPLE midpoint measurements")
    ax.loglog(
        steps,
        np.exp(intercept) * np.asarray(steps) ** slope,
        "--",
        color="#D55E00",
        label=f"log-log fit: order {slope:.3f}",
    )
    ax.set(
        xlabel="Step size (normalized time)",
        ylabel="Max absolute endpoint angle error (rad)",
        title="Analytic tokamak: midpoint vs refined Gauss2",
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(runner.output / "convergence.png", dpi=150)
    plt.close(fig)
    return [
        row(
            "Midpoint order of accuracy",
            "two-stage Gauss (order 4), 256 steps; stability checked at 128 steps",
            "1.8 <= fitted order <= 2.2; reference gap < 1% of smallest error",
            f"order={slope:.6g}; reference gap={reference_gap:.3g}; errors={errors}",
            1.8 <= slope <= 2.2 and reference_gap < 0.01 * min(errors),
        ),
        row(
            "Analytic tokamak invariants",
            "J_perp=p²(1-pitch²)/(1-r cos(theta)); P_phi=-r²/2+r⁴+p pitch(1+r cos(theta))",
            "max relative J_perp excursion <= 1e-10; max absolute native P_phi excursion <= 1e-10",
            f"J_perp={magnetic_moment_error:.6g}; P_phi={momentum_error:.6g}",
            magnetic_moment_error <= 1e-10 and momentum_error <= 1e-10,
        ),
        row(
            "Collisionless energy",
            "constant H0=(v/v0)² in the static analytic TEST field",
            "max relative excursion <= 1e-3 over normalized time 0.1",
            f"max relative excursion={energy:.6g} (all four midpoint traces)",
            energy <= 1e-3,
        ),
    ]


def uncertainty(runner: Runner, plt) -> list[str]:
    fixture = runner.output / "analytic-chartmap.nc"
    shutil.copyfile(REPO / "build/test/tests/test_boozer_chartmap.nc", fixture)
    rows, fractions = [], []
    for seed in (11, 22, 33):
        rng = np.random.default_rng(seed)
        starts = np.column_stack(
            [
                np.full(256, 0.85),
                rng.uniform(0, 2 * np.pi, 256),
                rng.uniform(0, np.pi, 256),
                np.ones(256),
                rng.uniform(-1, 1, 256),
            ]
        )
        _, _, outcomes = simple_run(
            runner,
            f"loss-seed-{seed}",
            dict(
                field_input="'../analytic-chartmap.nc'",
                coord_input="'../analytic-chartmap.nc'",
                isw_field_type=2,
                startmode=2,
                integmode=3,
                trace_time=1e-6,
                ntimstep=8,
                npoiper2=256,
                ntestpart=256,
                sbeg=0.85**2,
            ),
            starts,
        )
        if not np.array_equal(outcomes[:, 0], np.arange(1, 257)):
            raise ValueError("loss outcomes reordered or missing markers")
        lost = loss_mask(outcomes[:, 1], outcomes[:, 2], 1e-6)
        fractions.append(float(lost.mean()))
        for n in (32, 64, 128, 256):
            k = int(lost[:n].sum())
            lo, hi = wilson_interval(k, n)
            rows.append((seed, n, k, k / n, math.sqrt(k / n * (1 - k / n) / n), lo, hi))
    np.savetxt(
        runner.output / "loss-uncertainty.csv",
        rows,
        delimiter=",",
        header="seed,N,losses,fraction,binomial_SE,Wilson_low,Wilson_high",
        comments="",
    )
    lines = [
        "## Sampling uncertainty",
        "",
        "**Measured:** SIMPLE on the generated analytic Boozer chartmap; 256 independent starts per seed, "
        "rho=0.85, uniform theta in [0,2π), phi in [0,π), pitch in [-1,1], speed=1; "
        "trace_time=1e-6 s, midpoint, npoiper2=256. This is a diagnostic sampling law, not a reactor loss prediction.",
        "",
        "| Seed | N | Losses | Measured fraction | Binomial SE | 95% Wilson interval |",
        "|---|---|---|---|---|---|",
    ]
    for seed, n, k, p, se, lo, hi in rows:
        lines.append(
            f"| {seed} | {n} | {k} | {p:.6f} | {se:.6f} | [{lo:.6f}, {hi:.6f}] |"
        )
    lines += [
        "",
        f"**Measured seed spread at N=256:** sample SD={np.std(fractions, ddof=1):.6g}; "
        f"range={np.ptp(fractions):.6g}. Three seeds are descriptive, not a coverage test.",
        "",
        "**Calculated:** for iid Bernoulli losses, SE=sqrt(p(1-p)/N), at most 1/(2 sqrt(N)). "
        "Worst-case SE at N=100, 1000, 10000 is 0.05, 0.01581, 0.005. "
        f"[Wilson score formula]({NIST}). Prefixes within a seed are correlated; seeds are independent. "
        "Intervals quantify sampling only, conditional on the chosen field, time, step and start distribution. "
        "They exclude discretization bias, correlated/weighted production markers, equilibrium and model uncertainty. "
        "Unresolved/skipped markers fail the report and are never counted as confined or lost.",
        "",
    ]
    fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
    for seed, marker, color in zip(
        (11, 22, 33), ("o", "s", "^"), ("#0072B2", "#D55E00", "#009E73")
    ):
        data = np.array([r for r in rows if r[0] == seed])
        ax.errorbar(
            data[:, 1],
            data[:, 3],
            yerr=[data[:, 3] - data[:, 5], data[:, 6] - data[:, 3]],
            marker=marker,
            linestyle="--",
            capsize=3,
            color=color,
            label=f"seed {seed}",
        )
    ax.set(
        xscale="log",
        xlabel="Marker count N (nested prefixes)",
        ylabel="Measured loss fraction",
        title="SIMPLE losses with 95% Wilson intervals",
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(runner.output / "sampling.png", dpi=150)
    plt.close(fig)
    return lines


def provenance() -> list[str]:
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()

    cache = (REPO / "build/CMakeCache.txt").read_text()
    keys = ("CMAKE_BUILD_TYPE", "CMAKE_Fortran_COMPILER", "SIMPLE_DETERMINISTIC_FP")
    entries = [
        line
        for line in cache.splitlines()
        if any(line.startswith(k + ":") for k in keys)
    ]
    libneo = git("-C", "build/_deps/libneo-src", "rev-parse", "HEAD")
    compiler_match = re.search(r"^CMAKE_Fortran_COMPILER:[^=]+=(.+)$", cache, re.M)
    if compiler_match is None:
        raise ValueError("Fortran compiler missing from CMakeCache.txt")
    compiler = compiler_match[1]
    compiler_version = subprocess.check_output(
        [compiler, "--version"], text=True
    ).splitlines()[0]
    binary_hash = hashlib.sha256((REPO / "build/simple.x").read_bytes()).hexdigest()
    return [
        f"- UTC: {datetime.now(timezone.utc).isoformat()}",
        f"- Platform: {platform.platform()}; compiler: {compiler_version}",
        f"- simple.x SHA256: `{binary_hash}`",
        f"- SIMPLE: `{git('rev-parse', 'HEAD')}`; tracked changes present: {bool(git('diff', 'HEAD', '--stat'))}",
        f"- libneo: `{libneo}`",
        f"- Python: {sys.version.split()[0]}; NumPy: {np.__version__}",
        "- OMP_NUM_THREADS=1; OPENBLAS_NUM_THREADS=1",
        *[f"- `{entry}`" for entry in entries],
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("build/output/verification-report")
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=300,
        help="total subprocess budget in seconds (default 300)",
    )
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if not output.is_relative_to((REPO / "build").resolve()):
        parser.error("output directory must be inside this worktree's build/")
    if output.exists() and any(output.iterdir()):
        parser.error("use an empty output directory to avoid stale evidence")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be positive and finite")
    output.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lines = [
        "# SIMPLE numerical verification report",
        "",
        "Evidence for the requested EUROfusion criteria 3.1 and 3.4; this bounded suite is not a certification "
        "or experimental validation. Every numeric result below is measured unless marked calculated.",
        "",
        f"Command: `{shlex.join([sys.executable, *sys.argv])}`",
        "",
        "## Provenance",
        "",
    ]
    checks = []
    code = 1
    try:
        runner = Runner(output, args.timeout)
        checks.extend(consumer_checks(runner))
        lines += provenance()
        checks.extend(convergence(runner, plt))
        sampling = uncertainty(runner, plt)
        lines += [
            "",
            "## Checks",
            "",
            "| Check | Reference/oracle | Tolerance | Measured result | Status |",
            "|---|---|---|---|---|",
        ]
        for check in checks:
            # Escape pipes in scientific expressions so the Markdown remains a table.
            lines.append(
                "| "
                + " | ".join(str(value).replace("|", r"\|") for value in check.values())
                + " |"
            )
        lines += [
            "",
            "![Measured midpoint convergence](convergence.png)",
            "",
            "Endpoint error is max(|delta theta|, |delta phi|), using unwrapped radians in the analytic TEST field. "
            "Gauss radial and pitch outputs correspond to its last stage (src/orbit_symplectic.f90), "
            "so they are not used as endpoint reference values. The finite-time reference checks temporal order, not independent physics. "
            "All four fixed steps are fitted, with no point selection. Raw steps/errors are in convergence.csv.",
            "",
            *sampling,
            "![Measured loss sampling uncertainty](sampling.png)",
            "",
            "## Scope and evidence",
            "",
            "consumer-tests.log contains the registered tests and their full output. Each orbit/loss directory "
            "contains simple.in, start.dat, orbit_exit_code.dat, orbits.nc and run.log. "
            "The JSON table and CSV files preserve the reported measurements.",
            "",
            "Golden-record revision comparisons and the golden comparator sanity suite are **NOT RUN** here. "
            "Golden records require a second build and are covered separately by make test-golden / make test-regression. "
            "Long-time near-axis energy, SPECTRE, full-orbit, external equilibrium and inter-code/experimental "
            "comparisons remain outside this bounded report. See docs/verification.rst for the source inventory.",
            "",
        ]
        code = int(any(check["status"] == "FAIL" for check in checks))
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        lines += [
            "",
            "## FAIL: incomplete verification",
            "",
            str(exc),
            "",
            "A failed, missing, nonfinite or timed-out check is not a pass. Inspect the preserved logs.",
        ]
        print(str(exc), file=sys.stderr)
    (output / "measurements.json").write_text(
        json.dumps(checks, indent=2, allow_nan=False) + "\n"
    )
    (output / "verification-report.md").write_text("\n".join(lines) + "\n")
    print(f"{'FAIL' if code else 'PASS'}: {output / 'verification-report.md'}")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
