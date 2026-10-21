#!/usr/bin/env python3
"""Reject accidental bulk result artifacts before commit; never run tests."""

import os
import subprocess
import sys

MAX_FILES = 24
MAX_BYTES = 4 * 1024 * 1024
BINARY_SUFFIXES = (".npz", ".npy", ".gz", ".h5", ".hdf5", ".nc", ".zip", ".tar", ".zst")


def git(*args):
    return subprocess.check_output(("git", *args))


def is_result_artifact(path):
    lower = path.lower()
    return lower.endswith(BINARY_SUFFIXES) or (
        "/results/" in "/" + lower and lower.endswith((".csv", ".dat", ".jsonl"))
    )


def main():
    try:
        names = git("diff", "--cached", "--diff-filter=AM", "--name-only", "-z")
        artifacts = [os.fsdecode(p) for p in names.split(b"\0") if p and is_result_artifact(os.fsdecode(p))]
        sizes = [(p, int(git("cat-file", "-s", ":" + p))) for p in artifacts]
    except (subprocess.CalledProcessError, ValueError) as exc:
        print(f"staged artifact check failed: {exc}", file=sys.stderr)
        return 2

    total = sum(size for _, size in sizes)
    if len(sizes) <= MAX_FILES and total <= MAX_BYTES:
        return 0

    details = ", ".join(p for p, _ in sizes[:4])
    print(f"Bulk staged data: {len(sizes)} result files, {total / 1048576:.1f} MiB; "
          f"example paths: {details}", file=sys.stderr)
    print("Prefer an external, versioned data archive with checksums and small fixtures in Git.", file=sys.stderr)
    if os.environ.get("SIMPLE_ALLOW_BULK_ARTIFACTS") == "1":
        print("Explicit override SIMPLE_ALLOW_BULK_ARTIFACTS=1 accepted.", file=sys.stderr)
        return 0
    print("If this addition is scientifically necessary, retry with "
          "SIMPLE_ALLOW_BULK_ARTIFACTS=1 and justify it in the PR.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
