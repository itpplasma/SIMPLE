#!/usr/bin/env python3
"""The integration microstep dtaumin must depend on npoiper2 only, not trace_time.

Regression for the trace_time coupling bug: dtaumin used to be set to dtau/ntau,
which made the physical resolution depend on the total trace length. Two runs
differing only in trace_time then integrated at different resolution and their
confined-fraction curves were not comparable. dtaumin is now 2*pi*rbig/npoiper2.
"""

from __future__ import annotations

import math
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SIMPLE_EXE = REPO_ROOT / "build" / "simple.x"

NPOIPER2 = 256
NTIMSTEP = 100


def _write_simple_in(path: Path, vmec: Path, trace_time: float) -> None:
    (path / "simple.in").write_text(
        f"""&config
startmode = 1
ntestpart = 1
ntimstep = {NTIMSTEP}
trace_time = {trace_time}
npoiper2 = {NPOIPER2}
notrace_passing = 0
isw_field_type = -1
integmode = 3
netcdffile = '{vmec}'
deterministic = .true.
sbeg(1) = 0.35
num_surf = 1
/
"""
    )


def _run_dtaumin(tmp: Path, vmec: Path, trace_time: float) -> tuple[float, float]:
    _write_simple_in(tmp, vmec, trace_time)
    out = subprocess.run(
        [str(SIMPLE_EXE), "simple.in"], cwd=tmp, check=True,
        capture_output=True, text=True,
    ).stdout
    line = next(ln for ln in out.splitlines() if ln.strip().startswith("tau:"))
    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?", line.split("tau:")[1])
    dtaumin = float(nums[1].replace("d", "e").replace("D", "e"))
    last_t = np.loadtxt(tmp / "confined_fraction.dat")[-1, 0]
    return dtaumin, last_t


@pytest.mark.usefixtures("vmec_file")
def test_dtaumin_independent_of_trace_time(tmp_path: Path, vmec_file: str) -> None:
    if not SIMPLE_EXE.exists():
        pytest.skip("simple.x is not built; run CMake/Ninja build first")
    vmec = Path(vmec_file).resolve()

    short = tmp_path / "short"; short.mkdir()
    long = tmp_path / "long"; long.mkdir()
    dt_short, t_short = _run_dtaumin(short, vmec, 10.0)
    dt_long, t_long = _run_dtaumin(long, vmec, 100.0)

    # dtaumin identical across a 10x trace_time change (was ~2x apart before the fix)
    assert dt_short == pytest.approx(dt_long, rel=1e-12)
    # and equal to the resolution target 2*pi*rbig/npoiper2 (rbig=1 for the TEST field)
    assert dt_short == pytest.approx(2.0 * math.pi / NPOIPER2, rel=1e-12)
    # the macrostep schedule still spans each requested trace_time
    assert t_short == pytest.approx(10.0, rel=2e-2)
    assert t_long == pytest.approx(100.0, rel=2e-2)
