import subprocess
import sys
import tempfile
from pathlib import Path

CONFIG = """&config
nper = 2
npoiper = 16
ntimstep = 8
ntestpart = 4
trace_time = 1d-6
sbeg = 0.5d0
phibeg = 0.d0
thetabeg = 0.d0
contr_pp = -1.d0
facE_al = 100.d0
npoiper2 = 8
n_e = 2
n_d = 4
netcdffile = '{equilibrium}'
ns_s = 5
ns_tp = 5
multharm = 3
isw_field_type = 1
startmode = 1
integmode = 0
relerr = 1d-8
tcut = -1d0
deterministic = .true.
ran_seed = 12345
batch_size = 4
/
"""


def run(simple: Path, equilibrium: Path, directory: Path) -> None:
    config = directory / "simple.in"
    config.write_text(CONFIG.format(equilibrium=equilibrium), encoding="ascii")
    subprocess.run([simple, config], cwd=directory, check=True, timeout=30)


def read_table(path: Path) -> list[list[float]]:
    return [[float(value) for value in line.split()] for line in path.read_text().splitlines()]


def main() -> None:
    simple, wout, vmec_input = map(Path, sys.argv[1:4])
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        wout_run = root / "wout"
        direct_run = root / "direct"
        wout_run.mkdir()
        direct_run.mkdir()
        run(simple, wout, wout_run)
        run(simple, vmec_input, direct_run)

        wout_exit = read_table(wout_run / "orbit_exit_code.dat")
        direct_exit = read_table(direct_run / "orbit_exit_code.dat")
        assert [row[1] for row in direct_exit] == [row[1] for row in wout_exit]

        wout_loss = read_table(wout_run / "confined_fraction.dat")
        direct_loss = read_table(direct_run / "confined_fraction.dat")
        assert [row[1:3] for row in direct_loss] == [row[1:3] for row in wout_loss]

        wout_orbits = read_table(wout_run / "times_lost.dat")
        direct_orbits = read_table(direct_run / "times_lost.dat")
        assert len(direct_orbits) == len(wout_orbits)


if __name__ == "__main__":
    main()
