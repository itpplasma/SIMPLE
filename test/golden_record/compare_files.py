#!/usr/bin/env python3
"""Compare two numerical files for golden record tests.

Golden record tests primarily validate numerical reproducibility. Historically
we required bit-identical output, but refactors (e.g. spline kernels) can change
floating-point evaluation order while preserving physical correctness.

This comparator defaults to strict floating-point tolerances and can be tuned
via environment variables:
  - GOLDEN_RECORD_RTOL
  - GOLDEN_RECORD_ATOL
"""

from __future__ import annotations

import os
import sys

import numpy as np


def _env_float(name: str, default: float) -> float:
    value = os.environ.get(name, "")
    if value.strip() == "":
        return default
    return float(value)


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: compare_files.py <ref_file> <cur_file>")
        return 2

    ref_file = sys.argv[1]
    cur_file = sys.argv[2]

    rtol = _env_float("GOLDEN_RECORD_RTOL", 1.0e-7)
    atol = _env_float("GOLDEN_RECORD_ATOL", 1.0e-12)

    ref_data = np.loadtxt(ref_file)
    cur_data = np.loadtxt(cur_file)

    if ref_data.ndim == 1:
        ref_data = ref_data.reshape(-1, 1)
    if cur_data.ndim == 1:
        cur_data = cur_data.reshape(-1, 1)

    if ref_data.shape != cur_data.shape:
        print(
            f"Shape mismatch: {ref_file} has {ref_data.shape}, "
            f"{cur_file} has {cur_data.shape}"
        )
        return 1

    diff = np.abs(ref_data - cur_data)
    denom = np.maximum(np.abs(ref_data), 1.0e-300)
    rel = diff / denom
    ok = np.isclose(ref_data, cur_data, rtol=rtol, atol=atol)

    if bool(np.all(ok)):
        print(
            f"Files {ref_file} and {cur_file} match "
            f"(rtol={rtol:.1e}, atol={atol:.1e})."
        )
        return 0

    idx = np.argwhere(~ok)
    print(
        f"Found {idx.shape[0]} non-matching entries comparing {ref_file} and "
        f"{cur_file} (rtol={rtol:.1e}, atol={atol:.1e})."
    )
    print("First 5 mismatches (row, col, ref_val, cur_val, rel_diff):")
    for row, col in idx[:5]:
        print(
            f"  [{row},{col}]: ref={ref_data[row, col]:.16e}, "
            f"cur={cur_data[row, col]:.16e}, rel_diff={rel[row, col]:.2e}"
        )
    if idx.shape[0] > 5:
        print(f"  ... and {idx.shape[0] - 5} more")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
