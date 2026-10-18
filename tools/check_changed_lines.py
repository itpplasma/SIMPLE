#!/usr/bin/env python3
"""Check only added/modified lines, including staged and unstaged tracked edits."""

from __future__ import annotations

import argparse
import difflib
import json
import re
import subprocess
import sys
from pathlib import Path

HUNK_RE = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", re.MULTILINE)
FORTRAN_SUFFIXES = {".f90", ".f95", ".f03", ".f08", ".fpp"}


def run(
    command: list[str], *, check: bool = False, input: str | None = None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command, check=check, input=input, text=True, capture_output=True
    )


def changed_lines(base: str) -> dict[Path, set[int]]:
    """Diff the merge base against tracked working files, in source coordinates."""
    ancestor = run(["git", "merge-base", base, "HEAD"], check=True).stdout.strip()
    names = run(
        ["git", "diff", "--no-ext-diff", "--no-renames", "--name-only", "-z", ancestor],
        check=True,
    ).stdout
    files: dict[Path, set[int]] = {}
    for name in filter(None, names.split("\0")):
        path = Path(name)
        if not path.is_file():
            continue
        diff = run(
            [
                "git",
                "diff",
                "--no-ext-diff",
                "--no-textconv",
                "--no-renames",
                "--unified=0",
                ancestor,
                "--",
                name,
            ],
            check=True,
        ).stdout
        lines: set[int] = set()
        for match in HUNK_RE.finditer(diff):
            start = int(match.group(1))
            count = int(match.group(2) or "1")
            lines.update(range(start, start + count))
        if lines:
            files[path] = lines
    return files


def formatting_lines(source: str, formatted: str) -> set[int]:
    """Locate formatter edits in the ORIGINAL file, even after line wrapping."""
    original = source.splitlines(keepends=True)
    result: set[int] = set()
    matcher = difflib.SequenceMatcher(
        None, original, formatted.splitlines(keepends=True), autojunk=False
    )
    for tag, start, end, _, _ in matcher.get_opcodes():
        if tag in ("replace", "delete"):
            result.update(range(start + 1, end + 1))
        elif tag == "insert":
            # Attribute inserted whitespace to the following source line (or EOF).
            result.add(min(start + 1, len(original)))
    return result


def check_ruff(path: Path, lines: set[int]) -> list[str]:
    """Keep Ruff diagnostics whose primary source line was added/modified."""
    result = run(["ruff", "check", "--output-format=json", "--no-cache", str(path)])
    if result.returncode not in (0, 1):
        raise RuntimeError(f"ruff failed for {path}: {result.stderr.strip()}")
    diagnostics = json.loads(result.stdout)
    return [
        f"{path}:{item['location']['row']}:{item['location']['column']}: "
        f"{item['code']} {item['message']}"
        for item in diagnostics
        if item["location"]["row"] in lines
    ]


def check_formatter(command: list[str], path: Path, lines: set[int]) -> list[str]:
    """Compare formatted stdout without writing any file."""
    source = path.read_text()
    result = run([*command, str(path)], input=source)
    # fprettify 0.3.7 logs parse errors but can return success.
    if result.returncode or re.search(
        r"^(ERROR|CRITICAL):", result.stderr, re.MULTILINE
    ):
        raise RuntimeError(f"formatter failed for {path}: {result.stderr.strip()}")
    affected = formatting_lines(source, result.stdout) & lines
    return [f"{path}:{line}: formatting differs" for line in sorted(affected)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", default="origin/main", help="merge-base revision")
    args = parser.parse_args()
    failures: list[str] = []
    try:
        for path, lines in sorted(changed_lines(args.base).items()):
            if path.suffix.lower() in FORTRAN_SUFFIXES:
                failures.extend(
                    check_formatter(
                        ["fprettify", "--stdout", "-c", ".fprettify"], path, lines
                    )
                )
            elif path.parts[:1] in (("python",), ("tools",)) and path.suffix == ".py":
                failures.extend(check_ruff(path, lines))
                failures.extend(
                    check_formatter(
                        ["ruff", "format", "--no-cache", "--quiet", "--stdin-filename"],
                        path,
                        lines,
                    )
                )
    except (OSError, subprocess.CalledProcessError, ValueError, RuntimeError) as error:
        print(f"changed-line lint: {error}", file=sys.stderr)
        return 2

    if failures:
        print("Changed-line lint failed:\n" + "\n".join(failures))
        return 1
    print("Changed-line lint passed (only added/modified lines checked).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
