"""Compare native neutral-beam budgets for identical launches on one or two ranks."""

import math
import os
import subprocess
import sys
from pathlib import Path

MARKERS = 32
CONFIG = f"""&config
field_input = 'wout.nc'
coord_input = 'wout.nc'
isw_field_type = 1
integmode = 0
startmode = 8
nbi_file = 'beamlets.dat'
nbi_path_length = 20d0
nbi_torque_bins = 8
ntestpart = {MARKERS}
trace_time = 2d-5
ntimstep = 5
particle_energy_eV = 1d5
n_d = 2
n_e = 1
am1 = 2d0
am2 = 3d0
Z1 = 1d0
Z2 = 1d0
swcoll = .false.
contr_pp = -1d10
deterministic = .true.
ran_seed = 42
/
&profiles
Te_scale = 1d3
Ti1_scale = 1d3
ni1_scale = 3d19
ni2_scale = 0d0
/
"""
BEAMLETS = "T 2 1e5 1e6 0.5 0.3 0.2 6.0 -9.0 0.0 0.0 1.0 0.0\n"
TOTAL_KEYS = (
    "bulk_total",
    "field_total",
    "lost_angular_momentum_rate",
    "still_fast_angular_momentum_rate",
    "injected_ionized_angular_momentum_rate",
    "beam_target_neutron_rate_per_s",
)


def run_case(launcher, count_flag, simple, wout, case, ranks):
    case.mkdir(parents=True, exist_ok=True)
    wout_link = case / "wout.nc"
    wout_link.unlink(missing_ok=True)
    wout_link.symlink_to(wout)
    (case / "simple.in").write_text(CONFIG)
    (case / "beamlets.dat").write_text(BEAMLETS)
    command = [*launcher, count_flag, "1", simple, "simple.in"]
    if ranks == 2:
        # Separate seeds reproduce divergent sampled data, with identical orbit
        # dynamics once the root launches are broadcast.
        (case / "peer.in").write_text(
            CONFIG.replace("ran_seed = 42", "ran_seed = 2026")
        )
        command += [":", count_flag, "1", simple, "peer.in"]
    output = subprocess.run(
        command,
        cwd=case,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
        capture_output=True,
        text=True,
        timeout=120,
    )
    (case / "run.log").write_text(output.stdout + output.stderr)
    if output.returncode:
        raise RuntimeError(f"MPI run failed\n{output.stdout}\n{output.stderr}")
    if (
        ranks == 2
        and f"MPI ranks = 2, particles on this rank = {MARKERS // 2}"
        not in output.stdout
    ):
        raise AssertionError("MPI did not split the markers")


def read_budget(path):
    summary, bins = {}, []
    counts = ()
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            bins.append(tuple(map(float, line.split())))
        elif line.startswith("# markers begun "):
            counts = tuple(map(int, line.split()[-4:]))
        else:
            tokens = line[1:].split()
            if len(tokens) == 2 and tokens[0] in TOTAL_KEYS:
                summary[tokens[0]] = float(tokens[1])
    if len(counts) != 4 or len(bins) != 8:
        raise AssertionError(f"missing marker counts or radial bins in {path}")
    return counts, tuple(summary[key] for key in TOTAL_KEYS), bins


def check(case):
    counts, total, bins = read_budget(case / "nbi_torque.dat")
    if counts[0] != MARKERS or sum(counts[1:]) != MARKERS:
        raise AssertionError(f"MPI NBI marker counts are incomplete: {counts}")
    if not math.isclose(math.fsum(total[:4]), total[4], rel_tol=1e-12):
        raise AssertionError(f"MPI NBI budget does not close: {total[:5]}")
    if not total[-1] > 0:
        raise AssertionError("fixture must produce beam-target neutrons")
    print(f"{case.name}: counts {counts}, budget {total}")
    return counts, total, bins


def compare(rows_1, rows_2, label):
    if len(rows_1) != len(rows_2):
        raise AssertionError(f"{label} row count differs")
    for column, (values_1, values_2) in enumerate(
        zip(zip(*rows_1), zip(*rows_2), strict=True)
    ):
        scale = max(map(abs, values_1))
        for value_1, value_2 in zip(values_1, values_2, strict=True):
            if not math.isclose(value_1, value_2, rel_tol=1e-10, abs_tol=1e-12 * scale):
                raise AssertionError(
                    f"{label} column {column} differs: {value_1} vs {value_2}"
                )


def main():
    mpiexec, count_flag, simple, wout, directory = sys.argv[1:]
    version = subprocess.run([mpiexec, "--version"], capture_output=True, text=True)
    launcher = [mpiexec]
    if any(name in version.stdout + version.stderr for name in ("Open MPI", "OpenRTE")):
        launcher += ["--oversubscribe"]
    one, two = (Path(directory) / f"np{n}" for n in (1, 2))
    for case, n in ((one, 1), (two, 2)):
        run_case(launcher, count_flag, simple, Path(wout), case, n)
    counts_1, total_1, bins_1 = check(one)
    counts_2, total_2, bins_2 = check(two)
    if counts_1 != counts_2:
        raise AssertionError("MPI NBI marker outcomes differ")
    compare([total_1], [total_2], "NBI totals")
    compare(bins_1, bins_2, "NBI bins")
    for name in (
        "nbi_births.dat",
        "nbi_source.dat",
        "times_lost.dat",
        "orbit_exit_code.dat",
    ):
        if (one / name).read_bytes() != (two / name).read_bytes():
            raise AssertionError(f"MPI {name} differs")


if __name__ == "__main__":
    main()
