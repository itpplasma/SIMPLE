#!/usr/bin/env python3
"""Minimal examples of the SIMPLE Python API for orbit tracing and classification."""

from __future__ import annotations

from pathlib import Path

import pysimple


def trace_example(vmec_file: str | Path, n_particles: int = 32) -> None:
    """Trace particles and report confinement statistics."""

    # Initialize SIMPLE with VMEC file
    pysimple.init(vmec_file, deterministic=True, trace_time=5e-5, ntestpart=n_particles)

    # Sample particles on a flux surface
    particles = pysimple.sample_surface(n_particles, s=0.3)

    # Trace orbits in parallel
    results = pysimple.trace_parallel(particles, integrator="midpoint")

    # Report statistics
    confined = (results['loss_times'] >= 5e-5).sum()
    lost = (results['loss_times'] < 5e-5).sum()
    print(f"Traced {n_particles} particles: confined={confined}, lost={lost}")


def classify_example(vmec_file: str | Path, n_particles: int = 32) -> None:
    """Classify particle orbits as trapped/passing and regular/chaotic."""

    # Initialize with classification enabled (tcut > 0)
    pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, tcut=0.1, ntestpart=n_particles)

    # Sample particles
    particles = pysimple.sample_surface(n_particles, s=0.5)

    # Classify orbits
    results = pysimple.classify_parallel(particles)

    # Analyze classification
    n_passing = results['passing'].sum()
    n_trapped = (~results['passing']).sum()
    n_regular = (results['minkowski'] == 1).sum()
    n_chaotic = (results['minkowski'] == 2).sum()

    print(f"Classification of {n_particles} particles:")
    print(f"  Passing: {n_passing}, Trapped: {n_trapped}")
    print(f"  Regular: {n_regular}, Chaotic: {n_chaotic}")


def main() -> None:
    # Use test VMEC file
    repo_root = Path(__file__).resolve().parents[1]
    vmec_file = repo_root / "test" / "test_data" / "wout.nc"

    if not vmec_file.exists():
        print(f"VMEC file not found: {vmec_file}")
        print("Download with:")
        print("  wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O test/test_data/wout.nc")
        return

    print("=== Orbit Tracing Example ===")
    trace_example(vmec_file)

    print("\n=== Classification Example ===")
    classify_example(vmec_file)


if __name__ == "__main__":
    main()
