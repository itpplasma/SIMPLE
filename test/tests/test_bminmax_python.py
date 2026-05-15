#!/usr/bin/env python3
"""Backend-free Python tests for bminmax cache predicates."""

from __future__ import annotations

import sys
import importlib.util
from pathlib import Path


def main() -> int:
    helper = Path(__file__).resolve().parents[2] / "python" / "pysimple" / "_bminmax.py"
    spec = importlib.util.spec_from_file_location("pysimple_bminmax", helper)
    if spec is None or spec.loader is None:
        raise AssertionError(f"could not load {helper}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    needs_bminmax_cache = module.needs_bminmax_cache

    cases = [
        ((1, 0, False), False),
        ((0, 0, False), True),
        ((2, 0, False), True),
        ((0, 0, True), False),
        ((2, 0, True), True),
        ((0, 1, False), False),
        ((2, 1, False), True),
    ]

    for args, expected in cases:
        actual = needs_bminmax_cache(*args)
        if actual is not expected:
            raise AssertionError(
                f"needs_bminmax_cache{args} returned {actual}, expected {expected}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
