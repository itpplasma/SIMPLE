"""Classification by pysimple must match simple.x for the same simple.in."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

pysimple = pytest.importorskip("pysimple")

SIMPLE_X = Path(__file__).resolve().parents[2] / "build" / "simple.x"


def test_classification_matches_simple_x(tmp_path, vmec_file, monkeypatch):
    if not SIMPLE_X.exists():
        pytest.skip("simple.x is not built")
    (tmp_path / "simple.in").write_text(f"""&config
startmode = 1
ntestpart = 10
trace_time = 1d-3
npoiper2 = 64
tcut = 0.1
notrace_passing = 0
isw_field_type = 0
integmode = 3
netcdffile = '{Path(vmec_file).resolve()}'
deterministic = .true.
sbeg = 0.35
/
""")
    subprocess.run([str(SIMPLE_X), "simple.in"], cwd=tmp_path, check=True,
                   capture_output=True)
    # times_lost.dat: ipart, loss_time, trap_par, s_start, perp_inv, ...
    reference = np.loadtxt(tmp_path / "times_lost.dat")

    python_dir = tmp_path / "python"
    python_dir.mkdir()
    monkeypatch.chdir(python_dir)
    pysimple.init(config=tmp_path / "simple.in")
    result = pysimple.classify_parallel(pysimple.load_particles(tmp_path / "start.dat"))

    np.testing.assert_array_equal(result["passing"], reference[:, 2] < 0.0)
    for key, column in [("loss_times", 1), ("trap_parameter", 2),
                        ("perpendicular_invariant", 4)]:
        np.testing.assert_allclose(result[key], reference[:, column], rtol=1e-10,
                                   atol=1e-12, equal_nan=True, err_msg=key)
