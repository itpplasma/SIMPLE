"""Integration test: restart from partial times_lost.dat produces identical output."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

CONFIG_TEMPLATE = """\
&config
  isw_field_type = -1
  ntestpart = 32
  ntimstep = 50
  npoiper = 10
  trace_time = 1d-2
  deterministic = .true.
  integmode = 0
  sbeg = 0.3
  checkpoint_interval = 0.0
  {restart_line}
/
"""


def run_simple(simple_x: Path, work_dir: Path, timeout: int = 60) -> None:
    result = subprocess.run(
        [str(simple_x)],
        cwd=work_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=False,
    )
    (work_dir / "simple.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        raise AssertionError(
            f"simple.x failed (rc={result.returncode}):\n{result.stdout}"
        )


def read_dat(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if line.strip()]


def main() -> None:
    simple_x = Path(sys.argv[1]).resolve()
    work_root = Path(
        os.environ.get("CTEST_BINARY_DIRECTORY", Path.cwd())
    )

    # 1. Reference run
    ref_dir = Path(tempfile.mkdtemp(prefix="restart_ref_", dir=work_root))
    (ref_dir / "simple.in").write_text(
        CONFIG_TEMPLATE.format(restart_line=""), encoding="utf-8"
    )
    run_simple(simple_x, ref_dir)

    ref_tl = read_dat(ref_dir / "times_lost.dat")
    ref_cf = read_dat(ref_dir / "confined_fraction.dat")
    assert len(ref_tl) == 32, f"expected 32 particles, got {len(ref_tl)}"

    # 2. Build partial times_lost.dat: keep first 16, zero out rest
    rst_dir = Path(tempfile.mkdtemp(prefix="restart_rst_", dir=work_root))
    partial_lines = []
    for i, line in enumerate(ref_tl):
        cols = line.split()
        if i < 16:
            partial_lines.append(line)
        else:
            idx = cols[0]
            # Sentinel: times_lost=-1, zero trap_par/perp_inv/zend
            partial_lines.append(
                f"  {idx}  -1.00000000000000  0.00000000000000"
                f"  {cols[3]}"  # preserve zstart(1)
                f"  0.00000000000000"
                f"  0.00000000000000  0.00000000000000"
                f"  0.00000000000000  0.00000000000000"
                f"  0.00000000000000"
            )

    (rst_dir / "times_lost.dat").write_text(
        "\n".join(partial_lines) + "\n", encoding="utf-8"
    )
    (rst_dir / "simple.in").write_text(
        CONFIG_TEMPLATE.format(restart_line="restart = .true."),
        encoding="utf-8",
    )

    # 3. Restart run
    run_simple(simple_x, rst_dir)

    rst_tl = read_dat(rst_dir / "times_lost.dat")
    rst_cf = read_dat(rst_dir / "confined_fraction.dat")

    # 4. Compare times_lost.dat line by line
    assert len(rst_tl) == len(ref_tl), (
        f"row count mismatch: ref={len(ref_tl)} rst={len(rst_tl)}"
    )
    tl_mismatches = 0
    for i, (r, s) in enumerate(zip(ref_tl, rst_tl)):
        rv = [float(x) for x in r.split()]
        sv = [float(x) for x in s.split()]
        for j, (a, b) in enumerate(zip(rv, sv)):
            if abs(a - b) > 1e-10 * max(1.0, abs(a)):
                tl_mismatches += 1
                print(
                    f"times_lost mismatch line {i+1} col {j}: "
                    f"ref={a} rst={b}"
                )
    assert tl_mismatches == 0, f"{tl_mismatches} mismatches in times_lost.dat"

    # 5. Compare confined_fraction.dat
    assert len(rst_cf) == len(ref_cf), (
        f"confined_fraction row mismatch: ref={len(ref_cf)} rst={len(rst_cf)}"
    )
    cf_mismatches = 0
    for i, (r, s) in enumerate(zip(ref_cf, rst_cf)):
        rv = [float(x) for x in r.split()]
        sv = [float(x) for x in s.split()]
        for j, (a, b) in enumerate(zip(rv, sv)):
            if abs(a - b) > 1e-10 * max(1.0, abs(a)):
                cf_mismatches += 1
                print(
                    f"confined_fraction mismatch line {i+1} col {j}: "
                    f"ref={a} rst={b}"
                )
    assert cf_mismatches == 0, (
        f"{cf_mismatches} mismatches in confined_fraction.dat"
    )

    print("Restart integration test PASSED")


if __name__ == "__main__":
    main()
