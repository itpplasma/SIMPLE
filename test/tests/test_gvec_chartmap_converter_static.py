#!/usr/bin/env python3
"""Static checks for the GVEC chartmap converter CLI contract."""

from __future__ import annotations

import ast
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"FAIL: {message}")


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    script = repo / "tools" / "gvec_to_boozer_chartmap.py"
    tree = ast.parse(script.read_text(encoding="utf-8"))

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr != "add_argument":
            continue
        if not node.args or not isinstance(node.args[0], ast.Constant):
            continue
        if node.args[0].value != "--Bcov":
            continue
        for keyword in node.keywords:
            if keyword.arg == "default" and isinstance(keyword.value, ast.Constant):
                if keyword.value.value == "boozer-avg":
                    print("PASS: GVEC converter defaults to Boozer covariant averaging")
                    return
                fail(f"--Bcov default is {keyword.value.value!r}")
        fail("--Bcov has no default")

    fail("--Bcov argument not found")


if __name__ == "__main__":
    main()
