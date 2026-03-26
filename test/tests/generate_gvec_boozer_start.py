from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.optimize import root

import gvec


START_POINTS = np.array(
    [
        [0.35, 0.2, 0.1, 1.0, 0.2],
        [0.42, 1.0, 0.4, 1.0, -0.1],
        [0.58, 2.2, 0.7, 1.0, 0.4],
        [0.71, 4.1, 1.0, 1.0, -0.3],
    ],
    dtype=np.float64,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Boozer-coordinate SIMPLE start points from logical ones."
    )
    parser.add_argument("--param-file", required=True, type=Path)
    parser.add_argument("--state-file", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mnfactor", type=int, default=1)
    return parser.parse_args()


def angle_error(angle: float, target: float, period: float) -> float:
    return ((angle - target + 0.5 * period) % period) - 0.5 * period


def solve_boozer_point(
    state: gvec.State, s: float, theta: float, zeta: float, mnfactor: int
) -> tuple[float, float]:
    rho = np.array([np.sqrt(s)], dtype=np.float64)
    theta_period = 2.0 * np.pi
    zeta_period = 2.0 * np.pi / int(state.nfp)

    def residual(tz: np.ndarray) -> np.ndarray:
        ev = gvec.evaluate_sfl(
            state,
            "theta",
            "zeta",
            rho=rho,
            theta=np.array([tz[0]], dtype=np.float64),
            zeta=np.array([tz[1]], dtype=np.float64),
            sfl="boozer",
            MNfactor=mnfactor,
        )
        theta_eval = float(ev["theta"].values[0, 0, 0])
        zeta_eval = float(ev["zeta"].values[0, 0, 0])
        return np.array(
            [
                angle_error(theta_eval, theta, theta_period),
                angle_error(zeta_eval, zeta, zeta_period),
            ],
            dtype=np.float64,
        )

    theta_seeds = [theta, theta + 0.5 * np.pi, theta + np.pi, theta + 1.5 * np.pi]
    zeta_seeds = [
        zeta,
        zeta + 0.25 * zeta_period,
        zeta + 0.5 * zeta_period,
        zeta + 0.75 * zeta_period,
    ]
    best = None
    best_norm = np.inf

    for theta_seed in theta_seeds:
        for zeta_seed in zeta_seeds:
            guess = np.array(
                [theta_seed % theta_period, zeta_seed % zeta_period], dtype=np.float64
            )
            result = root(residual, guess, method="hybr", tol=1.0e-12)
            norm = np.linalg.norm(residual(result.x), ord=np.inf)
            if norm < best_norm:
                best = result.x.copy()
                best_norm = norm
            if result.success and norm < 1.0e-10:
                return float(result.x[0] % theta_period), float(result.x[1] % zeta_period)

    if best is None or best_norm >= 1.0e-8:
        raise RuntimeError(
            f"Failed to solve Boozer point for s={s}, theta={theta}, zeta={zeta}; residual={best_norm}"
        )
    return float(best[0] % theta_period), float(best[1] % zeta_period)


def main() -> None:
    args = parse_args()
    state = gvec.State(str(args.param_file), str(args.state_file))
    rows = []
    for s, theta, zeta, lam, vpar in START_POINTS:
        theta_b, zeta_b = solve_boozer_point(state, s, theta, zeta, args.mnfactor)
        rows.append([s, theta_b, zeta_b, lam, vpar])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output, np.asarray(rows, dtype=np.float64), fmt="%.16e")


if __name__ == "__main__":
    main()
