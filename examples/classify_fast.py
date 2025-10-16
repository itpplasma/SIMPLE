#!/usr/bin/env python3
"""Minimal example showing the fast Python classification API."""

from __future__ import annotations

from pathlib import Path

import numpy as np

import simple


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

    simple.load_vmec(vmec_file)

    batch = simple.ParticleBatch(16)
    batch.initialize_from_samplers(vmec_file, method="surface", s=0.4)

    result = simple.classify_fast(batch)

    labels = np.stack([result.j_parallel, result.topology, result.minkowski], axis=1)
    print("Classification codes (J_parallel, topology, Minkowski):")
    print(labels)


if __name__ == "__main__":
    main()
