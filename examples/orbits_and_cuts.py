#!/usr/bin/env python3
"""Trace one orbit and derive a simple toroidal-plane cut from the trajectory."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import pysimple
from simple_backend import get_can_sub


VMEC_FILE = Path(__file__).resolve().parents[1] / "test" / "test_data" / "wout.nc"
TRACE_TIME = 1.0e-2
INTEGRATOR = "midpoint"
INITIAL_POSITION = np.array([0.6, 0.0, 0.25 * np.pi, 1.0, -0.4], dtype=float)
TOROIDAL_CUT = 0.0
TOROIDAL_PERIOD = 2.0 * np.pi


def _valid_prefix_length(times: np.ndarray) -> int:
    if times.ndim != 1 or times.size == 0 or not np.isfinite(times[0]):
        return 0

    n_valid = 1
    for idx in range(1, times.size):
        if not np.isfinite(times[idx]) or times[idx] <= times[idx - 1]:
            break
        n_valid += 1
    return n_valid


def _find_plane_crossings(
    trajectory: np.ndarray,
    *,
    target: float,
    period: float,
) -> np.ndarray:
    if trajectory.shape[1] < 2:
        return np.empty((trajectory.shape[0], 0), dtype=np.float64)

    toroidal_turns = np.unwrap(trajectory[2, :]) - target
    normalized_turns = toroidal_turns / period
    cuts: list[np.ndarray] = []

    for idx in range(trajectory.shape[1] - 1):
        turn0 = normalized_turns[idx]
        turn1 = normalized_turns[idx + 1]
        if not np.isfinite(turn0) or not np.isfinite(turn1) or turn0 == turn1:
            continue

        low = min(turn0, turn1)
        high = max(turn0, turn1)
        first_boundary = int(np.floor(low)) + 1
        last_boundary = int(np.floor(high))

        for boundary in range(first_boundary, last_boundary + 1):
            alpha = (boundary - turn0) / (turn1 - turn0)
            if 0.0 <= alpha <= 1.0:
                cut = (1.0 - alpha) * trajectory[:, idx] + alpha * trajectory[:, idx + 1]
                cut[2] = target
                cuts.append(cut)

    if not cuts:
        return np.empty((trajectory.shape[0], 0), dtype=np.float64)
    return np.column_stack(cuts)


def _canonical_to_cyl(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rr = np.empty(points.shape[1], dtype=np.float64)
    pp = np.empty(points.shape[1], dtype=np.float64)
    zz = np.empty(points.shape[1], dtype=np.float64)

    for idx in range(points.shape[1]):
        s, theta_can, varphi_can = points[0, idx], points[1, idx], points[2, idx]
        theta_vmec, varphi_vmec = get_can_sub.can_to_vmec(
            float(s),
            float(theta_can),
            float(varphi_can),
        )
        rr[idx], zz[idx] = get_can_sub.vmec_to_cyl(
            float(s),
            float(theta_vmec),
            float(varphi_vmec),
        )
        rr[idx] /= 100.0
        zz[idx] /= 100.0
        pp[idx] = varphi_vmec

    return rr, pp, zz


def main() -> None:
    pysimple.init(
        VMEC_FILE,
        deterministic=True,
        trace_time=TRACE_TIME,
        ntestpart=1,
    )

    result = pysimple.trace_orbit(
        INITIAL_POSITION,
        integrator=INTEGRATOR,
        return_trajectory=True,
    )

    n_valid = _valid_prefix_length(result["times"])
    trajectory = result["trajectory"][:, :n_valid]
    if trajectory.shape[1] < 2:
        raise RuntimeError("Orbit trace did not produce enough trajectory samples.")

    valid_orbit = trajectory[:, (trajectory[0, :] > 0.0) & (trajectory[0, :] < 1.0)]
    if valid_orbit.shape[1] < 2:
        raise RuntimeError("Orbit left the valid flux domain before a cut could be computed.")

    cuts = _find_plane_crossings(
        valid_orbit,
        target=TOROIDAL_CUT,
        period=TOROIDAL_PERIOD,
    )
    if cuts.shape[1] == 0:
        raise RuntimeError("No toroidal-plane cuts were found; increase TRACE_TIME.")

    orbit_r, orbit_phi, orbit_z = _canonical_to_cyl(valid_orbit)
    cut_r, _, cut_z = _canonical_to_cyl(cuts)

    print(f"Traced orbit with {valid_orbit.shape[1]} retained samples")
    print(f"Computed {cuts.shape[1]} toroidal-plane cuts at phi = {TOROIDAL_CUT:.3f} rad")
    print(f"Final loss time estimate: {result['loss_time']:.3e} s")

    fig = plt.figure(figsize=(14, 6))
    ax_orbit = fig.add_subplot(121, projection="3d")
    ax_orbit.plot(
        orbit_r * np.cos(-orbit_phi),
        orbit_r * np.sin(-orbit_phi),
        orbit_z,
        linewidth=0.8,
        color="k",
    )
    ax_orbit.set_title("Orbit trajectory")
    ax_orbit.set_xlabel("x / m")
    ax_orbit.set_ylabel("y / m")
    ax_orbit.set_zlabel("z / m")

    ax_cut = fig.add_subplot(122)
    ax_cut.plot(cut_r, cut_z, ".", markersize=2, color="tab:red")
    ax_cut.set_title("Toroidal-plane cut")
    ax_cut.set_xlabel("R / m")
    ax_cut.set_ylabel("Z / m")
    ax_cut.axis("equal")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
