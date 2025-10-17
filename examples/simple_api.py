#!/usr/bin/env python3
"""Trace a small batch of particles with :class:`simple.SimpleSession`."""

from __future__ import annotations

from pathlib import Path

import simple


def run_trace_example(
    vmec_file: str | Path | None = None,
    *,
    n_particles: int = 32,
    surface: float = 0.3,
    tmax: float = 5e-5,
) -> simple.BatchResults:
    """Trace a batch of particles and return the raw :class:`BatchResults`."""

    vmec = vmec_file or simple.ensure_example_vmec()
    session = simple.SimpleSession(vmec)
    batch = session.sample_surface(n_particles, surface=surface)
    return session.trace(batch, tmax=tmax, integrator="symplectic_midpoint")


def main() -> None:
    results = run_trace_example()
    confined = results.confined_mask().sum()
    lost = results.lost_mask().sum()
    print(
        f"Traced {results.n_particles} particles for {results.tmax:.2e} s: "
        f"confined={confined}, lost={lost}"
    )


if __name__ == "__main__":
    main()
