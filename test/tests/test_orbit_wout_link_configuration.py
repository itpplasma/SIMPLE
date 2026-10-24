"""Check native configure-time link repair and idempotence without field data."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("--cmake", required=True)
parser.add_argument("--helper", type=Path, required=True)
args = parser.parse_args()

with tempfile.TemporaryDirectory(prefix="simple-orbit-link-") as directory:
    root = Path(directory)
    source = root / "source"
    source.mkdir()
    target = source / "fixture.wout"
    target.write_text("synthetic filesystem fixture\n")
    link = source / "examples/orbit_output_test/wout.nc"
    link.parent.mkdir(parents=True)
    (source / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.22)\n"
        "project(OrbitWoutLink LANGUAGES NONE)\n"
        'include("${LINK_HELPER}")\n'
        'simple_ensure_orbit_wout_link("${CMAKE_CURRENT_SOURCE_DIR}/fixture.wout" '
        '"${CMAKE_CURRENT_SOURCE_DIR}/examples/orbit_output_test/wout.nc")\n'
    )

    def configure():
        subprocess.run(
            [
                args.cmake,
                "-S",
                str(source),
                "-B",
                str(root / "build"),
                f"-DLINK_HELPER={args.helper.resolve()}",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert link.is_symlink() and os.readlink(link) == str(target)

    def identity():
        stat = link.lstat()
        return stat.st_ino, stat.st_mtime_ns, stat.st_ctime_ns

    configure()  # Missing link is created.
    first = identity()
    configure()
    assert identity() == first, "matching link was replaced during configuration"
    link.unlink()
    link.symlink_to(source / "wrong.wout")
    configure()  # A wrong link is repaired.
    repaired = identity()
    configure()
    assert identity() == repaired, "repaired link was replaced during configuration"
    link.unlink()
    configure()  # A removed link is recreated.
    print(
        json.dumps(
            {
                "missing_link": "repaired",
                "wrong_link": "repaired",
                "matching_link": "unchanged",
            }
        )
    )
