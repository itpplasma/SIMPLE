#!/usr/bin/env python3
import json
import pathlib
import shlex
import sys


def main():
    commands = json.loads(pathlib.Path(sys.argv[1]).read_text())
    dependency_source = pathlib.Path(sys.argv[2]).resolve()
    dependency_commands = []

    for entry in commands:
        source = pathlib.Path(entry["file"]).resolve()
        if dependency_source / "src" not in source.parents:
            continue
        if source.suffix.lower() != ".f90":
            continue
        dependency_commands.append(entry["command"])

    if not dependency_commands:
        raise SystemExit("no libneo Fortran compile command found")

    for command in dependency_commands:
        options = shlex.split(command)
        if "-ffast-math" in options:
            raise SystemExit("libneo deterministic build enables -ffast-math")
        if "-ffp-contract=fast" in options:
            raise SystemExit("libneo deterministic build enables fast contraction")
        if "-fno-fast-math" not in options:
            raise SystemExit("libneo deterministic build lacks -fno-fast-math")
        if "-ffp-contract=off" not in options:
            raise SystemExit("libneo deterministic build lacks -ffp-contract=off")


if __name__ == "__main__":
    main()
