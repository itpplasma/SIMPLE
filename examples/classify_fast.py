#!/usr/bin/env python3
"""Minimal example showing the fast classification API."""

from __future__ import annotations

from pathlib import Path

import simple


def run_fast_classification(
    vmec_file: str | None = None,
    *,
    n_particles: int = 16,
    classification_time: float = 0.015,
    include_minkowski: bool = False,
    output_dir: str | None = None,
) -> simple.ClassificationResult:
    """
    Classify a small batch of particles using the fast classifiers.

    Parameters
    ----------
    vmec_file:
        Optional VMEC equilibrium.  Falls back to :func:`simple.ensure_example_vmec`.
    n_particles:
        Number of particles sampled on the surface.
    classification_time:
        Duration of the fast classification trace.
    include_minkowski:
        Whether to run the Minkowski fractal-dimension classifier.  Disabling it keeps
        the classification focused on the fast J_parallel/topology heuristics.
    output_dir:
        Optional directory for legacy ``fort.*`` files (mirrors ``simple.x`` output).
    """

    vmec = vmec_file or simple.ensure_example_vmec()
    session = simple.SimpleSession(vmec)
    batch = session.sample_surface(n_particles, surface=0.4)

    output_path = None if output_dir is None else Path(output_dir)

    return session.classify_fast(
        batch,
        classification_time=classification_time,
        assume_passing_confined=True,
        legacy_files=output_path is not None,
        output_dir=output_path,
        include_minkowski=include_minkowski,
    )


def main() -> None:
    result = run_fast_classification()

    print("Classification counts (fast modes only):")
    for category, counts in result.counts().items():
        if category == "minkowski":
            continue
        print(f"  {category}: {counts}")


if __name__ == "__main__":
    main()
