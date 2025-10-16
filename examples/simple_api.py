#!/usr/bin/env python3
"""Example using the higher-level SIMPLE Python session API."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.append("../python")

import simple  # noqa: E402  pylint: disable=wrong-import-position


def main() -> None:
    vmec_file = simple.ensure_example_vmec()

    session = simple.SimpleSession(vmec_file)
    batch = session.sample_surface(100, surface=0.3)

    results = session.trace(batch, tmax=0.1, integrator="symplectic_midpoint")

    confined = results.confined()
    lost = results.lost()
    print(
        f"Results: {confined.shape[1]} confined, "
        f"{lost['loss_times'].size} lost out of {results.n_particles} particles"
    )

    loss_times = results.loss_times
    time_points = np.linspace(0.0, results.tmax, 100)
    confined_fraction = [
        np.sum(loss_times >= t) / loss_times.size for t in time_points
    ]

    plt.figure(figsize=(8, 6))
    plt.plot(time_points, confined_fraction)
    plt.xlabel("Time [normalized]")
    plt.ylabel("Confined fraction")
    plt.title("Particle confinement vs time")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
