"""Check native history, resolved denominators and preserved terminal events."""

import math
from pathlib import Path
import subprocess
import sys
import tempfile


def rows(path):
    return [[float(x) for x in line.split()] for line in path.read_text().splitlines()]


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    subprocess.run([sys.argv[1]], cwd=root, check=True, timeout=30)
    exits = rows(root / "orbit_exit_code.dat")
    assert [int(row[1]) for row in exits] == [104, 104, 0, 1]
    assert all(math.isnan(row[2]) for row in exits[:2])
    assert [row[2] for row in exits[2:]] == [2.0, 0.5]
    history = rows(root / "orbit_warning_holds.dat")
    assert history == [[1, 3, 2, 0, 2.0], [2, 2, 3, 1, 1.75]]
    confined = rows(root / "confined_fraction.dat")
    assert [int(row[3]) for row in confined] == [4, 3, 2]
    for row, expected in zip(confined, (1.0, 2 / 3, 0.5)):
        assert math.isclose(row[1], expected, abs_tol=1e-14)
    unresolved = rows(root / "unresolved_fraction.dat")
    assert [row[1] for row in unresolved] == [0.0, 0.25, 0.5]
    # Only the ordinary physical exit contributes to the physical-loss rate.
    assert rows(root / "avg_inverse_t_lost.dat") == [[2.0]]
    print("Python native history/accounting readback PASS")
