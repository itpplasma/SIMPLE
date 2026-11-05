"""Collisional markers stop as thermalized below thermal_stop_factor*T_i."""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

THERMALIZED = 5


def run(simple: Path, wout: Path, factor: float, work: Path) -> dict:
    work.mkdir()
    (work / "wout.nc").symlink_to(wout)
    (work / "simple.in").write_text(f"""&config
multharm = 3
contr_pp = -1e10
trace_time = 1d-4
sbeg = 0.3d0
ntestpart = 8
netcdffile = 'wout.nc'
isw_field_type = 2
deterministic = .True.
swcoll = .True.
thermal_stop_factor = {factor:.6e}
output_results_netcdf = .true.
/
""")
    with open(work / "simple.log", "w") as log:
        subprocess.run([str(simple)], cwd=work, stdout=log,
                       stderr=subprocess.STDOUT, timeout=600, check=True)
    with Dataset(work / "results.nc") as ds:
        return {k: np.array(ds[k][:]) for k in ("orbit_exit_code", "times_lost")}


def main() -> None:
    simple, wout = Path(sys.argv[1]), Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        # Far above any 3.5 MeV/T_i ratio: every marker stops after one step.
        stopped = run(simple, wout, 1e12, tmp / "stop")
        assert np.all(stopped["orbit_exit_code"] == THERMALIZED), stopped
        trace = np.nanmax(stopped["times_lost"])
        assert np.allclose(stopped["times_lost"], trace), "thermalized markers lost"
        free = run(simple, wout, 0.0, tmp / "free")
        assert not np.any(free["orbit_exit_code"] == THERMALIZED), free
    print("PASS: thermal stop")


if __name__ == "__main__":
    main()
