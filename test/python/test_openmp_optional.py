#!/usr/bin/env python3
"""Every use of omp_lib must be conditional so serial builds still compile.

ENABLE_OPENMP=OFF (and any compiler whose OpenMP runtime ships no omp_lib.mod)
drops the module, so an unconditional ``use omp_lib`` breaks the build outright.
Fortran's ``!$`` sentinel is the portable guard: the line is a comment unless
the compiler is in OpenMP mode.

The oracle here is the language rule, not this repository's current contents:
any file that references an omp_lib symbol outside a sentinel would fail to
compile without OpenMP, whether or not it exists today.
"""

import pathlib
import re
import sys

# Symbols that only exist when omp_lib is in scope.
OMP_SYMBOL = re.compile(r"\bomp_(?:get|set|in|init|destroy|test)_\w+\s*\(")
USE_OMP = re.compile(r"^\s*use\s+omp_lib\b", re.IGNORECASE)
SENTINEL = re.compile(r"^\s*!\$\s")


def offending_lines(path: pathlib.Path) -> list[tuple[int, str]]:
    bad = []
    for number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if SENTINEL.match(line):
            continue
        stripped = line.lstrip()
        if stripped.startswith("!"):
            continue
        if USE_OMP.match(line) or OMP_SYMBOL.search(line):
            bad.append((number, line.strip()))
    return bad


def main() -> None:
    root = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).parents[2]
    sources = [
        p
        for directory in ("src", "app", "test")
        for p in (root / directory).rglob("*.f90")
    ]
    if not sources:
        raise SystemExit(f"no Fortran sources found under {root}")

    failures = {}
    for path in sources:
        bad = offending_lines(path)
        if bad:
            failures[path.relative_to(root)] = bad

    if failures:
        for path, lines in sorted(failures.items()):
            for number, text in lines:
                print(f"{path}:{number}: unguarded OpenMP use: {text}")
        raise SystemExit(
            f"{sum(len(v) for v in failures.values())} unguarded omp_lib reference(s); "
            "prefix with the '!$' sentinel and provide a serial fallback"
        )
    print(f"checked {len(sources)} Fortran sources; all omp_lib use is guarded")


if __name__ == "__main__":
    main()
