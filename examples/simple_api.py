#!/usr/bin/env python3
"""
Example using the cleaned SIMPLE Python API.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.append("../python")

import simple  # noqa: E402  pylint: disable=wrong-import-position


def main() -> None:
    vmec_file = Path("wout.nc")
    if not vmec_file.exists():
        print(
            "Download VMEC file first:\n"
            "  wget https://github.com/hiddenSymmetries/simsopt/raw/master/"
            "tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc "
            "-O wout.nc"
        )
        return

    print("Loading VMEC equilibrium...")
    simple.load_vmec(vmec_file)

    print("Sampling particles on s=0.3 surface...")
    sampler = simple.SurfaceSampler(vmec_file)
    raw_particles = sampler.sample_surface_fieldline(100, s=0.3)
    batch = simple.ParticleBatch.from_fortran_arrays(raw_particles)

    print("Tracing orbits with symplectic midpoint integrator...")
    results = simple.trace_orbits(
        batch, tmax=0.1, integrator="symplectic_midpoint", verbose=False
    )

    confined = results.confined()
    lost = results.lost()
    print(
        f"Results: {confined.shape[1]} confined, "
        f"{lost['loss_times'].size} lost out of {results.n_particles}"
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
