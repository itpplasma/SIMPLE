#!/usr/bin/env python3
"""Trace orbits with macrostep output enabled."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import simple


def trace_macrostep_example(
    vmec_file: Optional[str] = None,
    *,
    batch: Optional[simple.ParticleBatch] = None,
    n_particles: int = 8,
    surface: float = 0.35,
    tmax: float = 1.0e-3,
    output_dir: str | Path = "orbit_macro_output",
) -> simple.BatchResults:
    """Trace a batch of particles and emit ``orbits.nc`` macrostep output."""

    vmec = vmec_file or simple.ensure_example_vmec()
    session = simple.SimpleSession(vmec)

    if batch is None:
        particles = session.sample_surface(n_particles, surface=surface)
    else:
        particles = batch

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    cwd = Path.cwd()
    try:
        os.chdir(output_path)
        with simple.macrostep_output(True):
            results = session.trace(particles, tmax=tmax)
    finally:
        os.chdir(cwd)

    return results


def main() -> None:
    results = trace_macrostep_example()
    print(
        f"Traced {results.n_particles} particles for {results.tmax:.3e}s. "
        "Macrostep output written to orbit_macro_output/orbits.nc"
    )


if __name__ == "__main__":
    main()
