#!/usr/bin/env python3
"""Minimal fast-classification script."""

from __future__ import annotations

from pathlib import Path

import pysimple


def classify_fast_example(vmec_file: str | Path, n_particles: int = 16) -> dict:
    """Run fast classification for a small batch."""
    pysimple.init(
        vmec_file,
        deterministic=True,
        trace_time=1e-4,
        tcut=0.1,
        ntestpart=n_particles,
    )

    particles = pysimple.sample_surface(n_particles, s=0.4)
    results = pysimple.classify_parallel(particles)

    n_passing = results['passing'].sum()
    n_trapped = (~results['passing']).sum()

    print(f"Fast classification of {n_particles} particles:")
    print(f"  Passing: {n_passing}")
    print(f"  Trapped: {n_trapped}")

    return results


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    vmec_file = repo_root / "test" / "test_data" / "wout.nc"

    if not vmec_file.exists():
        print(f"VMEC file not found: {vmec_file}")
        return

    classify_fast_example(vmec_file)


if __name__ == "__main__":
    main()
