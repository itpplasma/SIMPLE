#!/usr/bin/env python3
"""Select only the independently qualified Albert reporting context."""
import hashlib
import json
import shlex
from pathlib import Path
import subprocess
import sys


LIBNEO = "266f0a3c844315522d8b620fa4b96b69aedfeafc"
# Only this independently observed Release reference is qualified.
REFERENCE = "c75e5295dcf8af1382ec2985ce748da4cd079919"
INPUTS = {
    "simple.in": "c4cf23c11f13338ffe25ec553ac2a5738c21497b3d2be6f5350fbc0b70210080",
    "wout.nc": "36926d748099bc7abc138d82bdd4b447d2579f08f50189f27d441cf96c1b51e0",
    "coils.simple": "bb6e72654ed912c041ce97d71d1c736ddaf75cbc4a3ce31bad00190a860f9ea3",
}
TABLES = {
    "start.dat": "811a26a8989f7f49fb8ffcfae9fad496b24ed39946c7baf966c9ebe94bc6497d",
    "times_lost.dat": "1d4d6ac77bb28da16ffaae986225d18ddeb1299d3e0b61ca5f5d45efca3fa1e9",
    "orbit_exit_code.dat": "39e6ad03bb936fa9553fcf8fcffa365f611ecca1b22f9f85749c0066debd0711",
}


def git(root, *args):
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True, stderr=subprocess.DEVNULL
    ).strip()


def matches(root, expected):
    return all(hashlib.sha256((root / name).read_bytes()).hexdigest() == digest
               for name, digest in expected.items())


def qualified_dependency(build):
    cache = {}
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if "=" in line and ":" in line.split("=", 1)[0]:
            key, value = line.split("=", 1)
            cache[key.split(":", 1)[0]] = value
    if cache.get("CMAKE_BUILD_TYPE") != "Release":
        return False
    if cache.get("SIMPLE_DETERMINISTIC_FP") != "ON":
        return False
    if cache.get("LIBNEO_DETERMINISTIC_FP") != "ON":
        return False
    source = Path(cache["libneo_SOURCE_DIR"]).resolve()
    if git(source, "rev-parse", "HEAD") != LIBNEO:
        return False
    if git(source, "status", "--porcelain", "--untracked-files=no"):
        return False
    # Bind the resolved provider to the configured compilation, including
    # FetchContent source overrides. An environment pin alone is insufficient.
    commands = json.loads((build / "compile_commands.json").read_text())
    options = [shlex.split(entry["command"]) for entry in commands
               if source / "src" in Path(entry["file"]).resolve().parents
               and Path(entry["file"]).suffix.lower() == ".f90"]
    return bool(options) and all(
        "-fno-fast-math" in flags and "-ffp-contract=off" in flags
        and "-ffast-math" not in flags and "-ffp-contract=fast" not in flags
        for flags in options)


def main():
    reference, reference_build, current_build, reference_run, current_run = map(Path, sys.argv[1:])
    qualified = (
        git(reference, "rev-parse", "HEAD") == REFERENCE
        and not git(reference, "status", "--porcelain", "--untracked-files=no")
        and qualified_dependency(reference_build)
        and qualified_dependency(current_build)
        and matches(reference_run, INPUTS)
        and matches(current_run, INPUTS)
        and matches(reference_run, TABLES)
    )
    return 0 if qualified else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, KeyError, ValueError, subprocess.CalledProcessError) as error:
        print(f"Unqualified decoded reporting context: {error}", file=sys.stderr)
        raise SystemExit(1)
