#!/usr/bin/env python3
"""Trace orbits with macrostep output enabled."""

from __future__ import annotations

from pathlib import Path

import pysimple


def trace_macrostep_example(
    vmec_file: str | Path,
    n_particles: int = 8,
) -> dict:
    """Trace particles with single orbit trajectory output."""
    pysimple.init(
        vmec_file,
        deterministic=True,
        trace_time=1e-4,
        ntestpart=1,
    )

    particles = pysimple.sample_surface(1, s=0.35)
    particle = particles[:, 0]

    result = pysimple.trace_orbit(particle, integrator="midpoint", return_trajectory=True)

    n_timesteps = result['times'].shape[0]
    print(f"Traced 1 particle with {n_timesteps} timesteps")
    print(f"Loss time: {result['loss_time']:.3e}s")

    return result


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    vmec_file = repo_root / "test" / "test_data" / "wout.nc"

    if not vmec_file.exists():
        print(f"VMEC file not found: {vmec_file}")
        return

    trace_macrostep_example(vmec_file)


if __name__ == "__main__":
    main()
