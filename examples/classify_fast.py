#!/usr/bin/env python3
"""Minimal fast-classification script (Minkowski disabled)."""

from __future__ import annotations

import simple


def classify_fast_example(vmec_file: str | None = None) -> simple.ClassificationResult:
    """Run the fast topological and J_parallel classifiers for a small batch."""
    vmec = vmec_file or simple.ensure_example_vmec()
    session = simple.SimpleSession(vmec)
    batch = session.sample_surface(16, surface=0.4)
    return session.classify_fast(batch, classification_time=0.015)


def main() -> None:
    result = classify_fast_example()
    print("Fast classifier counts:")
    for category, counts in result.counts().items():
        print(f"  {category}: {counts}")


if __name__ == "__main__":
    main()
