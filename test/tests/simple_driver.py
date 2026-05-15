"""Shared helpers for SIMPLE executable driver tests."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def run_simple_case(
    simple_x: Path,
    wout: Path,
    name: str,
    config: str,
    *,
    timeout: int = 90,
) -> Path:
    work_root = Path(os.environ.get("CTEST_BINARY_DIRECTORY", Path.cwd()))
    case_dir = Path(tempfile.mkdtemp(prefix=f"{name}_", dir=work_root))
    (case_dir / "simple.in").write_text(config, encoding="utf-8")
    link_or_copy(wout, case_dir / "wout.nc")

    result = subprocess.run(
        [str(simple_x)],
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=False,
    )
    (case_dir / "simple.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        raise AssertionError(
            f"{name}: simple.x failed with {result.returncode}\n{result.stdout}"
        )
    return case_dir


def link_or_copy(source: Path, target: Path) -> None:
    try:
        target.symlink_to(source)
    except OSError:
        shutil.copy2(source, target)


def log_text(case_dir: Path) -> str:
    return (case_dir / "simple.log").read_text(encoding="utf-8")


def assert_log_order(case_dir: Path, before: str, after: str) -> None:
    log = log_text(case_dir)
    before_idx = log.find(before)
    after_idx = log.find(after)
    if before_idx < 0 or after_idx < 0:
        raise AssertionError(f"{case_dir.name}: missing log markers\n{log}")
    if before_idx > after_idx:
        raise AssertionError(
            f"{case_dir.name}: wrong log order for {before!r} and {after!r}"
        )


def assert_file_absent(case_dir: Path, filename: str) -> None:
    if (case_dir / filename).exists():
        raise AssertionError(f"{case_dir.name}: {filename} should not be written")


def assert_row_count(case_dir: Path, filename: str, expected: int) -> None:
    path = case_dir / filename
    if not path.exists():
        raise AssertionError(f"{case_dir.name}: expected {filename}")
    rows = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(rows) != expected:
        raise AssertionError(
            f"{case_dir.name}: expected {expected} {filename} rows, got {len(rows)}"
        )
