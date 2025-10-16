#!/usr/bin/env python3
"""Minimal example showing the fast classification API."""

from __future__ import annotations

from pathlib import Path

import numpy as np

import simple


def main() -> None:
    vmec = simple.ensure_example_vmec()
    session = simple.SimpleSession(vmec)

    batch = session.sample_surface(16, surface=0.4)

    result = session.classify_fast(
        batch,
        classification_time=0.1,
        assume_passing_confined=True,
        legacy_files=False,
    )

    labels = np.vstack([result.j_parallel, result.topology, result.minkowski]).T
    print("Classification codes (J_parallel, topology, Minkowski):")
    print(labels)

    print("Counts:")
    for category, counts in result.counts().items():
        print(f"  {category}: {counts}")

    # Write legacy fort.* outputs for comparison
    output_dir = Path("classification_output")
    session.classify_fast(
        batch,
        classification_time=0.1,
        assume_passing_confined=True,
        legacy_files=True,
        output_dir=output_dir,
    )
    print(f"Legacy fort.* outputs written to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
