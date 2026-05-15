#!/usr/bin/env python3
"""Backend-free Python tests for bminmax cache predicates."""

from __future__ import annotations

import sys
import importlib.util
from pathlib import Path


def load_needs_bminmax_cache():
    helper = Path(__file__).resolve().parents[2] / "python" / "pysimple" / "_bminmax.py"
    spec = importlib.util.spec_from_file_location("pysimple_bminmax", helper)
    if spec is None or spec.loader is None:
        raise AssertionError(f"could not load {helper}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.needs_bminmax_cache


def iter_cases(path: Path):
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        num_surf, ntcut, class_plot, expected = line.split()
        yield (
            (int(num_surf), int(ntcut), class_plot == "T"),
            expected == "T",
        )


def main() -> int:
    needs_bminmax_cache = load_needs_bminmax_cache()
    cases_path = Path(sys.argv[1])

    for args, expected in iter_cases(cases_path):
        actual = bool(needs_bminmax_cache(*args))
        if actual is not expected:
            raise AssertionError(
                f"needs_bminmax_cache{args} returned {actual}, expected {expected}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
