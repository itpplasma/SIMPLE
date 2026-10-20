"""Reproducible IDS orbit comparisons; all datasets/results live outside Git.

Run from the repository root after installing pysimple[imas] and ASCOT5:
  python examples/imas_benchmark.py --case circular --output /tmp/ids-benchmark
Public equilibria are fetched by immutable URL and verified before use.
"""

import argparse
import hashlib
import json
from pathlib import Path
import time
import urllib.request
import re
import subprocess
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import numpy as np

# Published example equilibria: (immutable URL, SHA256, COCOS of the file).
# Downloaded on first use into --cache; nothing is stored in Git.
SOURCES = {
    # DIII-D shot 145419 at 2100 ms (EFIT), OMAS sample data, MIT licence.
    # Conventional single-null tokamak of ASDEX Upgrade size; EFIT stores |q|.
    "diii-d": (
        "https://raw.githubusercontent.com/gafusion/omas/"
        "f11358ecb710be8db84ab08e16e4bf468f52ff8b/omas/samples/g145419.02100",
        "087aefddacac4337d54347e1e73085ef3b21c254176885726841a4521174f81f",
        1,
    ),
    # ITER hybrid scenario (Citrin et al 2010 Nucl. Fusion 50 115007; CRONOS,
    # CHEASE), TORAX third-party data under LICENSE_PINT (MIT).
    "iter-hybrid": (
        "https://raw.githubusercontent.com/google-deepmind/torax/"
        "beab14be873d272f94c97aeb741ff093ccab8746/torax/data/third_party/geo/"
        "iterhybrid_cocos17.eqdsk",
        "02d948b51525156d883e5afd7085b7c5b68304341be49360d88405d3151ac367",
        17,
    ),
    # SPARC equilibrium-library case 1300 (0.2 MA scan), SPARCPublic, MIT.
    "sparc": (
        "https://raw.githubusercontent.com/cfs-energy/SPARCPublic/"
        "5b913b8216d05346b22e2e1e87c2ff100761613f/EQlibrary/sparc_1300.EQDSK_PsiOver2pi_negIp_negBt_negFpol",
        "e7054c18fd7684d22c8eb52c2c53bc2c305803bb7aee9090d83573d45733edc2",
        3,
    ),
}


def revision(checkout):
    """Record a checkout's revision, or explicit unknowns for installed packages."""

    def git(*args):
        return subprocess.check_output(
            ["git", "-C", str(checkout), *args], text=True, stderr=subprocess.DEVNULL
        ).strip()

    try:
        if (
            Path(git("rev-parse", "--show-toplevel")).resolve()
            != Path(checkout).resolve()
        ):
            return {"commit": None, "dirty": None}
        return {
            "commit": git("rev-parse", "HEAD"),
            "dirty": bool(git("status", "--porcelain")),
        }
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "dirty": None}


def fetch(case, cache):
    url, checksum, _ = SOURCES[case]
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{case}.geqdsk"
    data = (
        path.read_bytes()
        if path.exists()
        else urllib.request.urlopen(url, timeout=60).read()
    )
    if hashlib.sha256(data).hexdigest() != checksum:
        raise ValueError(f"Checksum mismatch: {path}")
    if not path.exists():
        temporary = path.with_suffix(f".download-{os.getpid()}")
        temporary.write_bytes(data)
        temporary.replace(path)
    return path


def equilibrium_from_geqdsk(path, cocos, *, version="4.1.1"):
    """Benchmark glue for explicitly identified right-handed COCOS conventions.

    The original R,phi,Z axes are right handed; no sign is guessed from an IDS.
    Both DD3 (COCOS 11) and DD4 (COCOS 17) have sigma_rho_theta_phi = +1, so
    sign(q) = sign(Ip B0); EFIT-style files that store |q| are thereby covered.
    Profiles and flux are converted independently of either code's IDS reader.
    """
    from libneo.eqdsk_base import read_eqdsk
    from scipy.interpolate import CubicSpline
    import imas
    from pysimple.imas._common import initialize

    sigma_bp = {1: 1, 3: -1, 5: 1, 7: -1, 11: 1, 13: -1, 15: 1, 17: -1}
    if cocos not in sigma_bp:
        raise ValueError("This benchmark supports right-handed COCOS 1/3/5/7/1x")
    d = read_eqdsk(path)
    eq = imas.IDSFactory(version=version).new("equilibrium")
    initialize(eq, [0.0])
    eq.vacuum_toroidal_field.r0 = d["R0"]
    eq.vacuum_toroidal_field.b0 = [d["Btor_at_R0"]]
    eq.time_slice.resize(1)
    t = eq.time_slice[0]
    t.time = 0.0
    # Target sigma_Bp is +1 for DD3 (COCOS 11) and -1 for DD4 (COCOS 17).
    scale = (2 * np.pi if cocos < 10 else 1) * sigma_bp[cocos]
    if not version.startswith("3."):
        scale = -scale
    if np.ptp(np.sign(d["qprof"])) != 0:
        raise ValueError(f"q changes sign in {path}")
    p1 = t.profiles_1d
    p1.psi = scale * np.linspace(d["PsiaxisVs"], d["PsiedgeVs"], d["nrgr"])
    p1.f, p1.pressure, p1.q = (
        d["fprof"],
        d["ptotprof"],
        np.sign(d["Ip"] * d["Btor_at_R0"]) * np.abs(d["qprof"]),
    )
    from scipy.integrate import cumulative_trapezoid

    toroidal = cumulative_trapezoid(np.asarray(p1.q), np.asarray(p1.psi), initial=0)
    p1.rho_tor_norm = np.sqrt(toroidal / toroidal[-1])
    t.profiles_2d.resize(1)
    p2 = t.profiles_2d[0]
    p2.grid_type.index = 1
    p2.grid.dim1, p2.grid.dim2 = d["R"], d["Z"]
    p2.psi = scale * d["PsiVs"].T
    norm = np.clip(
        (d["PsiVs"].T - d["PsiaxisVs"]) / (d["PsiedgeVs"] - d["PsiaxisVs"]), 0, 1
    )
    bphi = (
        CubicSpline(np.linspace(0, 1, len(d["fprof"])), d["fprof"])(norm)
        / d["R"][:, None]
    )
    setattr(p2, "b_field_tor" if version.startswith("3.") else "b_field_phi", bphi)
    q = t.global_quantities
    q.psi_axis, q.psi_boundary = scale * d["PsiaxisVs"], scale * d["PsiedgeVs"]
    q.magnetic_axis.r, q.magnetic_axis.z = d["Rpsi0"], d["Zpsi0"]
    q.ip = d["Ip"]
    boundary = d["Lcfs"]
    if not np.array_equal(boundary[0], boundary[-1]):
        boundary = np.vstack((boundary, boundary[0]))
    t.boundary.outline.r, t.boundary.outline.z = boundary.T
    eq.validate()
    return eq


def circular_equilibrium(*, aspect=10.0, sign=1, version="4.1.1"):
    """Exact circular flux surfaces; psi=C*r_minor^2 and F=R0*B0.

    Used as a magnetic test field, without claiming scalar-pressure MHD balance.
    Large aspect ratio is an explicit parameter, rather than a renamed tokamak.
    """
    import imas
    from pysimple.imas._common import initialize

    eq = imas.IDSFactory(version=version).new("equilibrium")
    initialize(eq, [0.0])
    eq.vacuum_toroidal_field.r0 = aspect
    eq.vacuum_toroidal_field.b0 = [2.0]
    eq.time_slice.resize(1)
    t = eq.time_slice[0]
    t.time = 0.0
    scale = (-1 if version.startswith("3.") else 1) * 2 * np.pi
    r = np.linspace(aspect - 1.3, aspect + 1.3, 129)
    z = np.linspace(-1.3, 1.3, 129)
    c = sign * 0.5
    t.profiles_2d.resize(1)
    p2 = t.profiles_2d[0]
    p2.grid_type.index = 1
    p2.grid.dim1, p2.grid.dim2 = r, z
    p2.psi = scale * c * ((r[:, None] - aspect) ** 2 + z[None, :] ** 2)
    setattr(
        p2,
        "b_field_tor" if version.startswith("3.") else "b_field_phi",
        2 * aspect / r[:, None] * np.ones((1, len(z))),
    )
    s = np.linspace(0, 1, len(r))
    p1 = t.profiles_1d
    p1.psi = scale * c * s
    p1.f, p1.pressure = np.full(len(s), 2 * aspect), np.zeros(len(s))
    p1.q = -2 * aspect / (2 * c * np.sqrt(aspect**2 - s))
    p1.rho_tor_norm = np.sqrt(
        (aspect - np.sqrt(aspect**2 - s)) / (aspect - np.sqrt(aspect**2 - 1))
    )
    q = t.global_quantities
    q.psi_axis, q.psi_boundary = 0, scale * c
    q.magnetic_axis.r, q.magnetic_axis.z = aspect, 0
    q.ip = -4 * np.pi * c / np.sqrt(aspect**2 - 1) / (4 * np.pi * 1e-7)
    theta = np.linspace(0, 2 * np.pi, 129)
    t.boundary.outline.r, t.boundary.outline.z = aspect + np.cos(theta), np.sin(theta)
    eq.validate()
    return eq


def field_consistency(equilibrium, path):
    """Check chart geometry and |B| against the original rectangular IDS grid."""
    import netCDF4
    from scipy.interpolate import RectBivariateSpline

    eq = equilibrium.time_slice[0]
    p = eq.profiles_2d[0]
    psi = RectBivariateSpline(p.grid.dim1, p.grid.dim2, p.psi)
    toroidal = RectBivariateSpline(p.grid.dim1, p.grid.dim2, p.b_field_phi)
    with netCDF4.Dataset(path) as dataset:
        x, y, z = [np.asarray(dataset[k][0]).T * 0.01 for k in ["x", "y", "z"]]
        radius = np.hypot(x, y)
        flux = psi.ev(radius, z)
        rho = np.asarray(dataset["rho"][:])
        interior = (rho >= 0.1) & (rho <= 0.9)
        spread = np.ptp(flux[interior], axis=1) / abs(
            float(eq.global_quantities.psi_boundary)
            - float(eq.global_quantities.psi_axis)
        )
        expected = np.sqrt(
            toroidal.ev(radius, z) ** 2
            + (psi.ev(radius, z, dx=1) ** 2 + psi.ev(radius, z, dy=1) ** 2)
            / (2 * np.pi * radius) ** 2
        )
        actual = np.asarray(dataset["Bmod"][0]).T * 1e-4
        return {
            "rho_range": [0.1, 0.9],
            "maximum_relative_flux_surface_spread": float(np.max(spread)),
            "maximum_relative_Bmod_error": float(
                np.max(np.abs(actual[interior] / expected[interior] - 1))
            ),
        }


def bulk_plasma():
    import imas
    from pysimple.imas._common import initialize

    cp = imas.IDSFactory(version="4.1.1").new("core_profiles")
    initialize(cp, [0.0])
    cp.profiles_1d.resize(1)
    p = cp.profiles_1d[0]
    p.time = 0.0
    p.grid.rho_tor_norm = np.linspace(0, 1, 101)
    p.electrons.density, p.electrons.temperature = np.full(101, 1e20), np.full(101, 1e4)
    p.ion.resize(2)
    for ion, mass in zip(p.ion, [2, 3]):
        ion.element.resize(1)
        ion.element[0].a, ion.element[0].z_n, ion.z_ion = mass, 1, 1
        ion.density_thermal, ion.temperature = np.full(101, 5e19), np.full(101, 1e4)
    cp.validate()
    return cp


def ascot_run(equilibrium, markers, path, duration, step, plasma=None):
    """Use ASCOT's actual IMAS readers and compiled guiding-centre integrator."""
    from a5py import Ascot
    from a5py.ascot5io.options import Opt
    import unyt

    a = Ascot(str(path), create=True)
    a.data.create_input("imas b2ds", equilibrium_ids=equilibrium)
    m = a.data.create_input(
        "imas marker", distribution_sources_ids=markers, dryrun=True
    )
    for name in [
        "E_TC",
        "plasma_1D",
        "N0_1D",
        "wall_2D",
        "Boozer",
        "MHD_STAT",
        "asigma_loc",
    ]:
        if name != "plasma_1D" or plasma is None:
            a.data.create_input(name)
    if plasma is not None:
        bfield = a.data.bfield.active.read()
        profiles = a.data.create_input(
            "imas plasma",
            equilibrium_ids=equilibrium,
            core_profiles_ids=plasma,
            psi0=bfield["psi0"],
            psi1=bfield["psi1"],
            dryrun=True,
        )
        # ASCOT applies collisions after each orbit step and before the RHOMAX
        # check; plasma_1D raises an input-evaluation error for rho >= its last
        # grid point, aborting markers that step across rho = 1 instead of
        # ending them with RHOMAX. Hold the edge values constant outside the
        # IDS grid, as ASCOT's slowing-down tutorial does with rho in [0, 2].
        profiles["rho"] = np.append(profiles["rho"], 2.0)
        profiles["nrho"] += 1
        for key in ["vtor", "edensity", "etemperature", "itemperature", "idensity"]:
            values = np.asarray(profiles[key])
            profiles[key] = np.concatenate((values, values[-1:]), axis=0)
        a.data.create_input("plasma_1D", **profiles)
    a.simulation_initinputs()
    options = Opt.get_default()
    if not np.isclose(duration / step, round(duration / step), rtol=0, atol=1e-6):
        raise ValueError("duration must be an integer multiple of ASCOT's fixed step")
    options.update(
        SIM_MODE=2,
        ENABLE_ORBIT_FOLLOWING=1,
        ENABLE_ADAPTIVE=0,
        FIXEDSTEP_USE_USERDEFINED=1,
        FIXEDSTEP_USERDEFINED=step,
        ADAPTIVE_TOL_ORBIT=1e-11,
        ADAPTIVE_MAX_DRHO=0.05,
        ADAPTIVE_MAX_DPHI=1,
        ENDCOND_SIMTIMELIM=1,
        # ASCOT checks time > limit after a completed step.
        ENDCOND_LIM_SIMTIME=duration - step / 2,
        ENDCOND_RHOLIM=1,
        ENDCOND_MAX_RHO=1.0,
        ENDCOND_MIN_RHO=0.0,
    )
    options["ENABLE_COULOMB_COLLISIONS"] = int(plasma is not None)
    a.simulation_initoptions(**options)
    # Match SIMPLE's integer-multiple-of-proton-mass, nonrelativistic ion model.
    # ASCOT's isotope number and SIMPLE's n_d both remain the IDS mass number.
    mass_number = float(markers.source[0].species.ion.element[0].a)
    m["mass"] = np.full(m["n"], mass_number * 1.6726e-27) * unyt.kg
    # Gyrocentre IDS markers have no gyrophase; ASCOT's simulation API asks
    # for one even for a gyroaveraged calculation, where it has no effect.
    m["zeta"] = np.zeros(m["n"]) * unyt.rad
    a.simulation_initmarkers(**m)
    start = time.perf_counter()
    result = a.simulation_run(printsummary=False)
    elapsed = time.perf_counter() - start
    r, phi, z, energy, pitch, condition, elapsed_time, rho, *error = result.getstate(
        "r",
        "phi",
        "z",
        "ekin",
        "pitch",
        "endcond",
        "time",
        "rho",
        "errormsg",
        "errorline",
        "errormod",
        state="end",
        mode="gc",
    )
    out = np.column_stack(
        (
            r.to_value("m"),
            phi.to_value("rad"),
            z.to_value("m"),
            energy.to_value("eV"),
            pitch.to_value("1"),
        )
    )
    codes = np.asarray(condition).astype(int).copy()
    times = elapsed_time.to_value("s").copy()
    native = {
        "rho": np.asarray(rho).astype(float).copy(),
        "error": np.column_stack([np.asarray(e).astype(int) for e in error]),
    }
    a.simulation_free(inputs=True, markers=True, diagnostics=True)
    return out, codes, times, elapsed, native


def gorilla_run(checkout, field, positions, directory, duration, resolution, species):
    """Compare GORILLA's native polynomial pusher on the same chartmap.

    One marker per native plotting invocation keeps trajectory records and
    invariants unambiguous. Input/output and logs remain in the run directory.
    GORILLA does not export final pitch through this plotting executable.
    """
    import pysimple

    checkout = checkout.expanduser().resolve()
    executable = checkout / "BUILD/test_gorilla_main.x"
    references = pysimple.cylindrical_to_reference(positions[:, :3].T).T
    final, times, energy_errors, momentum_errors = [], [], [], []
    begin = time.perf_counter()
    for i, (position, reference) in enumerate(zip(positions, references)):
        work = directory / f"gorilla-{i:03d}"
        work.mkdir()
        configs = {
            "gorilla.inp": {
                "coord_system": 2,
                "ispecies": species,
                "ipusher": 2,
                "poly_order": 4,
                "i_time_tracing_option": 2,
                "eps_Phi": 0.0,
            },
            "tetra_grid.inp": {
                "grid_kind": 6,
                "netcdf_filename": str(field),
                "n1": resolution,
                "n2": 2 * resolution,
                "n3": 12,
                "sfc_s_min": 0.01,
            },
            "gorilla_plot.inp": {
                "i_orbit_options": 1,
                "total_orbit_time": duration,
                "energy_eV_start": position[3],
                "boole_full_orbit": True,
                "n_skip_full_orbit": 1,
                "boole_poincare_phi_0": False,
                "boole_poincare_vpar_0": False,
                "boole_e_tot": True,
                "boole_p_phi": True,
                "boole_J_par": False,
                "filename_orbit_start_pos_sthetaphi": "start.dat",
                "filename_full_orbit_sthetaphi": "orbit.dat",
                "filename_e_tot": "energy.dat",
                "filename_p_phi": "momentum.dat",
            },
        }
        for filename, overrides in configs.items():
            text = (checkout / "INPUT" / filename).read_text()
            for name, value in overrides.items():
                encoded = (
                    (".true." if value else ".false.")
                    if isinstance(value, bool)
                    else ("'" + value + "'" if isinstance(value, str) else str(value))
                )
                pattern = rf"(?im)^(\s*{name}\s*=)[^\n]*"
                text, count = re.subn(
                    pattern, lambda m: m[1] + " " + encoded + ",", text
                )
                if count != 1:
                    raise ValueError(
                        f"Expected one GORILLA setting {name}, got {count}"
                    )
            (work / filename).write_text(text)
        np.savetxt(
            work / "start.dat",
            [[reference[0] ** 2, reference[1], reference[2], position[4]]],
        )
        with (work / "run.log").open("w") as log:
            subprocess.run(
                [str(executable)],
                cwd=work,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                timeout=120,
            )
        orbit = np.atleast_2d(np.loadtxt(work / "orbit.dat"))[-1].copy()
        orbit[0] = np.sqrt(orbit[0])
        final.append(pysimple.reference_to_cylindrical(orbit))
        energy = np.atleast_2d(np.loadtxt(work / "energy.dat"))
        momentum = np.atleast_2d(np.loadtxt(work / "momentum.dat"))
        times.append(float(energy[-1, 0]))
        energy_errors.append(float((energy[-1, 1] - energy[0, 1]) / energy[0, 1]))
        momentum_errors.append(
            float((momentum[-1, 1] - momentum[0, 1]) / momentum[0, 1])
        )
    return {
        "final_cylindrical": np.asarray(final).tolist(),
        "end_times_s": times,
        "relative_energy_drift": energy_errors,
        "relative_toroidal_momentum_drift": momentum_errors,
        "seconds_including_mesh_initialization": time.perf_counter() - begin,
    }


def cylinder(directory, duration, step, energy, species):
    """Uniform axial field: exact guiding-centre translation and gyro-motion.

    ASCOT uses its analytic field and IDS markers. An infinite straight
    cylinder has no closed toroidal surfaces for SIMPLE's IDS converter.
    """
    from a5py import Ascot
    from a5py.ascot5io.options import Opt
    from pysimple.imas import Markers, write_markers
    import unyt

    mass, charge = (4, 2) if species == "alpha" else (2, 1)
    source = write_markers(
        [
            Markers(
                np.array([[1, 0.3, 0, energy, 0.6]]),
                np.ones(1),
                mass,
                charge,
                charge,
                rate=True,
            )
        ]
    )
    report = {
        "case": "cylinder",
        "field": "uniform Bx=2 T; zero E (straight cylinder axis X)",
        "duration_s": duration,
        "step_s": step,
        "simple": "unsupported: the IDS converter requires closed toroidal flux surfaces",
    }
    for mode, name in [(2, "guiding_centre"), (1, "full_orbit")]:
        a = Ascot(str(directory / f"{name}.h5"), create=True)
        a.data.create_input(
            "B_TC",
            bxyz=np.array([[2.0], [0.0], [0.0]]),
            jacobian=np.zeros((3, 3)),
            rhoval=0.5,
        )
        for inp in [
            "E_TC",
            "plasma_1D",
            "N0_1D",
            "wall_2D",
            "Boozer",
            "MHD_STAT",
            "asigma_loc",
        ]:
            a.data.create_input(inp)
        m = a.data.create_input(
            "imas marker", distribution_sources_ids=source, dryrun=True
        )
        m["mass"] = np.full(m["n"], mass * 1.6726e-27) * unyt.kg
        a.simulation_initinputs()
        options = Opt.get_default()
        options.update(
            SIM_MODE=mode,
            ENABLE_ORBIT_FOLLOWING=1,
            ENABLE_ADAPTIVE=0,
            FIXEDSTEP_USE_USERDEFINED=1,
            FIXEDSTEP_USERDEFINED=step,
            ENDCOND_SIMTIMELIM=1,
            ENDCOND_LIM_SIMTIME=duration - step / 2,
        )
        a.simulation_initoptions(**options)
        a.simulation_initmarkers(**m)
        result = a.simulation_run(printsummary=False)
        coordinates = "gc" if mode == 2 else "prt"
        r, phi, z, vr, vp, vz = result.getstate(
            "r", "phi", "z", "vr", "vphi", "vz", state="ini", mode=coordinates
        )
        r, phi, z = (
            float(r.to_value("m")[0]),
            float(phi.to_value("rad")[0]),
            float(z.to_value("m")[0]),
        )
        vr, vp, vz = [float(v.to_value("m/s")[0]) for v in [vr, vp, vz]]
        x0, y0 = r * np.cos(phi), r * np.sin(phi)
        vx, vy = (
            vr * np.cos(phi) - vp * np.sin(phi),
            vr * np.sin(phi) + vp * np.cos(phi),
        )
        if mode == 2:
            vpar = float(
                result.getstate("vpar", state="ini", mode="gc").to_value("m/s")[0]
            )
            expected = [x0 + vpar * duration, y0, z]
        else:
            gamma = 1 / np.sqrt(1 - (vx * vx + vy * vy + vz * vz) / 299792458.0**2)
            omega = charge * 1.602176634e-19 * 2 / (gamma * mass * 1.6726e-27)
            angle = omega * duration
            expected = [
                x0 + vx * duration,
                y0 + (vy * np.sin(angle) + vz * (1 - np.cos(angle))) / omega,
                z + (vz * np.sin(angle) - vy * (1 - np.cos(angle))) / omega,
            ]
        r, phi, z, t, e = result.getstate(
            "r", "phi", "z", "time", "ekin", state="end", mode=coordinates
        )
        r, phi, z = (
            float(r.to_value("m")[0]),
            float(phi.to_value("rad")[0]),
            float(z.to_value("m")[0]),
        )
        actual = np.array([r * np.cos(phi), r * np.sin(phi), z])
        report[name] = {
            "position_error_m": float(np.linalg.norm(actual - expected)),
            "end_time_s": float(t.to_value("s")[0]),
            "relative_energy_drift": float(e.to_value("eV")[0] / energy - 1),
        }
        a.simulation_free(inputs=True, markers=True, diagnostics=True)
    (directory / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(report, indent=2))


def suite(args, output):
    """Run a matrix in isolated processes because both solvers hold global state."""
    jobs = []
    if args.suite in ("verification", "all"):
        for case in ["circular", "large-aspect", *SOURCES]:
            for resolution, step in [(48, 1e-8), (96, 1e-8), (96, 5e-9)]:
                options = [
                    "--case",
                    case,
                    "--resolution",
                    str(resolution),
                    "--step",
                    str(step),
                ]
                if args.gorilla and step == 1e-8:
                    options += ["--gorilla", str(args.gorilla)]
                jobs.append((f"{case}-r{resolution}-dt{step:g}", options))
        for step in [1e-10, 5e-11]:
            jobs.append(
                (
                    f"cylinder-dt{step:g}",
                    ["--case", "cylinder", "--duration", "1e-6", "--step", str(step)],
                )
            )
        for case in ["circular", "large-aspect", *SOURCES]:
            jobs.append(
                (
                    f"alpha-{case}",
                    [
                        "--case",
                        case,
                        "--species",
                        "alpha",
                        "--energy",
                        "3500000",
                        "--dense",
                        "--duration",
                        "1e-4",
                    ],
                )
            )
    if args.suite in ("collisions", "all"):
        for case in ["circular", *SOURCES]:
            for step in [1e-8, 5e-9]:
                jobs.append(
                    (
                        f"collisions-{case}-dt{step:g}",
                        [
                            "--case",
                            case,
                            "--collisions",
                            "--species",
                            "alpha",
                            "--energy",
                            "3500000",
                            "--dense",
                            "--copies",
                            str(max(20, args.copies)),
                            "--duration",
                            "1e-4",
                            "--step",
                            str(step),
                        ],
                    )
                )
    if args.suite in ("regimes", "all"):
        # Regime-aware collisionless comparison: self-divergence horizons,
        # pointwise checks for regular orbits and paired cohort statistics.
        for case in ["circular", "large-aspect", *SOURCES]:
            for species, energy, extra in [
                ("deuterium", "10000", []),
                ("alpha", "3500000", ["--stat-markers", "512"]),
            ]:
                jobs.append(
                    (
                        f"regimes-{species}-{case}",
                        [
                            "--case",
                            case,
                            "--regimes",
                            "--species",
                            species,
                            "--energy",
                            energy,
                            "--dense",
                            "--resolution",
                            "96",
                            "--duration",
                            "1e-4",
                            "--step",
                            "5e-9",
                            *extra,
                        ],
                    )
                )
    summary = {"suite": args.suite, "planned": len(jobs), "runs": []}

    def execute(job):
        name, options = job
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--output",
            str(output / name),
            "--cache",
            str(args.cache),
            *options,
        ]
        try:
            with (output / f"{name}.log").open("w") as log:
                result = subprocess.run(
                    command, stdout=log, stderr=subprocess.STDOUT, timeout=1800
                )
            return {
                "name": name,
                "exit_code": result.returncode,
                "report": str(output / name / "report.json"),
            }
        except subprocess.TimeoutExpired:
            return {"name": name, "exit_code": 124, "error": "Run exceeded 1800 s"}

    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        pending = [executor.submit(execute, job) for job in jobs]
        for future in as_completed(pending):
            result = future.result()
            summary["runs"].append(result)
            (output / "suite.json").write_text(json.dumps(summary, indent=2) + "\n")
            print(
                f"{len(summary['runs'])}/{len(jobs)} {result['name']}: exit {result['exit_code']}",
                flush=True,
            )
    if any(r["exit_code"] != 0 for r in summary["runs"]):
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=["cylinder", "circular", "large-aspect", *SOURCES],
        default="circular",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--cache", type=Path, default=Path.home() / ".cache/simple/imas-benchmarks"
    )
    parser.add_argument("--duration", type=float, default=1e-5)
    parser.add_argument("--step", type=float, default=1e-8)
    parser.add_argument("--resolution", type=int, default=48)
    parser.add_argument(
        "--suite", choices=["verification", "collisions", "regimes", "all"]
    )
    parser.add_argument(
        "--jobs", type=int, default=2, help="Concurrent isolated benchmark processes"
    )
    parser.add_argument(
        "--gorilla", type=Path, help="Built GORILLA checkout (optional third code)"
    )
    parser.add_argument(
        "--species", choices=["deuterium", "alpha"], default="deuterium"
    )
    parser.add_argument(
        "--energy", type=float, default=1e4, help="Initial kinetic energy [eV]"
    )
    parser.add_argument(
        "--dense",
        action="store_true",
        help="54 markers spanning trapped/passing and edge orbits",
    )
    parser.add_argument(
        "--collisions",
        action="store_true",
        help="Flat 10 keV DT plasma through core_profiles IDS; alpha markers only",
    )
    parser.add_argument(
        "--copies",
        type=int,
        default=1,
        help="Independent copies for collision statistics; input weight unchanged",
    )
    parser.add_argument(
        "--regimes",
        action="store_true",
        help="Collisionless regime-aware comparison (imas_trajectory_regimes.py)",
    )
    parser.add_argument(
        "--stat-markers",
        type=int,
        default=0,
        help="Random cohort size for --regimes loss statistics",
    )
    parser.add_argument(
        "--perturbation",
        type=float,
        default=1e-6,
        help="Start displacement in R [m] for --regimes self-divergence",
    )
    parser.add_argument(
        "--horizon-tolerance",
        type=float,
        default=1e-3,
        help="Self-divergence [m] that ends pointwise comparison (--regimes)",
    )
    args = parser.parse_args()
    if args.regimes and (args.collisions or args.case == "cylinder"):
        parser.error("--regimes compares collisionless toroidal orbits")
    if args.copies < 1 or (
        args.collisions and (args.species != "alpha" or args.case == "cylinder")
    ):
        parser.error("copies must be positive; native collisions require alpha markers")
    output = args.output.expanduser().resolve()
    repository = Path(__file__).resolve().parents[1]
    if output == repository or repository in output.parents:
        parser.error("Choose an output directory outside the repository")
    args.cache = args.cache.expanduser().resolve()
    if args.cache == repository or repository in args.cache.parents:
        parser.error("Choose a download cache outside the repository")
    for name in ("step", "duration", "energy", "perturbation", "horizon_tolerance"):
        value = getattr(args, name)
        if not np.isfinite(value) or value <= 0:
            parser.error(f"{name} must be positive and finite")
    if args.resolution < 8 or args.stat_markers < 0 or args.jobs < 1:
        parser.error("resolution must be >=8, stat-markers >=0 and jobs positive")
    if not np.isclose(
        args.duration / args.step, round(args.duration / args.step), rtol=0, atol=1e-6
    ):
        parser.error(
            "Require positive duration/step and an integer number of ASCOT steps"
        )
    output.mkdir(parents=True, exist_ok=False)
    if args.suite:
        return suite(args, output)
    if args.case == "cylinder":
        return cylinder(output, args.duration, args.step, args.energy, args.species)
    import imas
    import pysimple
    import a5py
    import libneo
    from pysimple.imas import Markers, equilibrium_to_chartmap, run, write_markers

    revisions = {
        "benchmark": revision(repository),
        "SIMPLE": revision(Path(pysimple.__file__).resolve().parents[2]),
        "ASCOT5": revision(Path(a5py.__file__).resolve().parents[1]),
        "libneo_converter": revision(Path(libneo.__file__).resolve().parents[2]),
    }
    if args.case in SOURCES:
        eq = equilibrium_from_geqdsk(
            fetch(args.case, args.cache), SOURCES[args.case][2]
        )
    else:
        eq = circular_equilibrium(aspect=100.0 if args.case == "large-aspect" else 10.0)
    t = eq.time_slice[0]
    r0, z0 = (
        float(t.global_quantities.magnetic_axis.r),
        float(t.global_quantities.magnetic_axis.z),
    )
    minor = max(t.boundary.outline.r) - r0
    radii = [0.1, 0.3, 0.5, 0.7, 0.9, 0.97] if args.dense else [0.3, 0.6, 0.8]
    pitches = np.linspace(-0.9, 0.9, 9) if args.dense else [-0.8, -0.2, 0.2, 0.8]
    positions = np.array(
        [[r0 + f * minor, 0.3, z0, args.energy, p] for f in radii for p in pitches]
    )
    positions = np.repeat(positions, args.copies, axis=0)
    mass, charge = (4, 2) if args.species == "alpha" else (2, 1)
    markers = write_markers(
        [
            Markers(
                positions,
                np.ones(len(positions)) / args.copies,
                mass,
                charge,
                charge,
                rate=True,
            )
        ]
    )
    plasma = bulk_plasma() if args.collisions else None
    with imas.DBEntry(str(output / "input.nc"), "w", dd_version="4.1.1") as entry:
        entry.put(eq)
        entry.put(markers)
        if plasma is not None:
            entry.put(plasma)
    n = args.resolution
    field = equilibrium_to_chartmap(
        eq,
        output / "field.nc",
        nstep=720,
        nlabel=2 * n,
        ntheta_int=2 * n,
        nsurfmax=2000,
        nsurf=4 * n,
        nrho=n,
        ntheta=2 * n,
    )
    if args.regimes:
        import imas_trajectory_regimes as regimes

        report = regimes.compare(
            eq,
            field,
            positions,
            mass,
            charge,
            duration=args.duration,
            step=args.step,
            output=output,
            perturbation=args.perturbation,
            horizon_tolerance=args.horizon_tolerance,
            stat_markers=args.stat_markers,
        )
        report.update(
            case=args.case,
            species=args.species,
            initial_energy_eV=args.energy,
            resolution=n,
            source=SOURCES.get(args.case),
            field_consistency_with_input_IDS=field_consistency(eq, field),
            revisions=revisions,
        )
        regimes.write(report, output)
        (output / "report.json").write_text(
            json.dumps(
                {k: v for k, v in report.items() if k != "markers"},
                indent=2,
                allow_nan=False,
            )
            + "\n"
        )
        print(json.dumps(report["summary"], indent=2))
        return
    start = time.perf_counter()
    actor_input = (
        {
            "equilibrium": eq,
            "core_profiles": plasma,
            "conversion": dict(
                nstep=720,
                nlabel=2 * n,
                ntheta_int=2 * n,
                nsurfmax=2000,
                nsurf=4 * n,
                nrho=n,
                ntheta=2 * n,
            ),
        }
        if plasma is not None
        else {"field": field}
    )
    result = run(
        markers,
        **actor_input,
        duration=args.duration,
        source_interval=1.0,
        workdir=output,
        integmode=pysimple.RK45,
        ntimstep=10000,
        deterministic=True,
        ran_seed=42,
    )
    simple_elapsed = time.perf_counter() - start
    simple = np.asarray(result.distribution[0].markers[1].positions)[:, :5]
    ascot, codes, times, ascot_elapsed, native = ascot_run(
        eq, markers, output / "ascot.h5", args.duration, args.step, plasma
    )
    from a5py.ascot5io.state import State

    simple_codes = np.asarray(result.distribution[0].markers[1].positions)[:, 6].astype(
        int
    )
    # Each marker must end with exactly one resolved condition: TLIM (confined
    # for the whole duration) or RHOMAX (lost). Both loss summaries use this.
    ascot_finished = np.array([State.endcond_check(int(c), "TLIM") for c in codes])
    ascot_lost = np.array([State.endcond_check(int(c), "RHOMAX") for c in codes])
    ascot_aborted = np.array(
        [State.endcond_check(int(c), "ABORTED") for c in codes]
    ) | (native["error"][:, 0] != 0)
    ascot_unresolved = ascot_aborted | (ascot_finished == ascot_lost)
    ascot_states = {}
    for c in codes:
        name = State.endcond_tostring(int(c)) or f"UNKNOWN({int(c)})"
        ascot_states[name] = ascot_states.get(name, 0) + 1
    matched = (
        (simple_codes == 0)
        & ascot_finished
        & ~ascot_unresolved
        & np.isclose(times, args.duration, rtol=1e-10, atol=1e-14)
    )
    jointly_completed = matched.copy()
    if args.collisions:
        matched[:] = False

    def cart(p):
        return np.column_stack(
            (p[:, 0] * np.cos(p[:, 1]), p[:, 0] * np.sin(p[:, 1]), p[:, 2])
        )

    report = {
        "case": args.case,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "revisions": revisions,
        "dd_version": "4.1.1",
        "duration_s": args.duration,
        "resolution": n,
        "ascot_fixedstep_s": args.step,
        "species": args.species,
        "collisions": args.collisions,
        "marker_copies": args.copies,
        "initial_energy_eV": args.energy,
        "simple_integrator": "RK45",
        "simple_ntimstep": 10000,
        "ascot_integrator": "fixed-step GC RK4",
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS", "runtime default"),
        "source": SOURCES.get(args.case),
        "field_consistency_with_input_IDS": field_consistency(eq, field),
        "markers": len(positions),
        "simple_seconds_including_field_initialization": simple_elapsed,
        "ascot_seconds_tracing": ascot_elapsed,
        "simple_exit_conditions": simple_codes.tolist(),
        "simple_lost_fraction": float(np.mean(simple_codes != 0)),
        "ascot_lost_fraction": float(np.mean(ascot_lost))
        if not ascot_unresolved.any()
        else None,
        "ascot_loss_definition": "RHOMAX; every other marker must be TLIM",
        "ascot_end_state_counts": ascot_states,
        "ascot_unresolved_markers": int(np.sum(ascot_unresolved)),
        "ascot_native_errors": [
            {
                "marker": int(i),
                "errormsg": int(native["error"][i, 0]),
                "errorline": int(native["error"][i, 1]),
                "errormod": int(native["error"][i, 2]),
                "rho": float(native["rho"][i]),
                "time_s": float(times[i]),
            }
            for i in np.flatnonzero(native["error"].any(axis=1))
        ],
        "jointly_completed": jointly_completed.tolist(),
        "trajectory_comparison": "disabled for independent collision streams"
        if args.collisions
        else "matched final times only",
        "position_error_m": [
            float(e) if ok else None
            for e, ok in zip(
                np.linalg.norm(cart(simple) - cart(ascot), axis=1), matched
            )
        ],
        "relative_energy_error": [
            float(e) if ok else None
            for e, ok in zip((simple[:, 3] - ascot[:, 3]) / positions[:, 3], matched)
        ],
        "pitch_error": [
            float(e) if ok else None
            for e, ok in zip(simple[:, 4] - ascot[:, 4], matched)
        ],
        "ascot_end_conditions": codes.tolist(),
        "ascot_end_times_s": times.tolist(),
    }
    if args.gorilla and not args.collisions:
        report["gorilla"] = gorilla_run(
            args.gorilla,
            field,
            positions,
            output,
            args.duration,
            max(12, n // 2),
            3 if mass == 4 else 2,
        )
        report["revisions"]["GORILLA"] = revision(args.gorilla)
    if args.gorilla and args.collisions:
        report["gorilla"] = "Not compared: GORILLA's core pusher is collisionless"

    def moments(end, completed, lost):
        survival = float(np.mean(completed))

        def error(values):
            if args.copies < 2:
                return None
            groups = values.reshape(-1, args.copies)
            return float(
                np.sqrt(
                    np.sum(np.var(groups, axis=1, ddof=1))
                    / (len(groups) ** 2 * args.copies)
                )
            )

        mean_pitch = (
            float(np.sum(end[:, 4] * completed) / np.sum(completed))
            if survival
            else None
        )
        return {
            "lost_fraction": float(np.mean(lost)),
            "lost_fraction_standard_error": error(lost.astype(float)),
            "remaining_energy_eV_per_initial_particle": float(
                np.mean(end[:, 3] * completed)
            ),
            "remaining_energy_standard_error_eV": error(end[:, 3] * completed),
            "confined_mean_pitch": mean_pitch,
            "confined_mean_squared_pitch": float(np.mean(end[completed, 4] ** 2))
            if survival
            else None,
            "confined_mean_energy_eV": float(np.mean(end[completed, 3]))
            if survival
            else None,
            "confined_energy_standard_error_eV": error(
                completed * (end[:, 3] - np.mean(end[completed, 3]))
            )
            / survival
            if survival and args.copies > 1
            else None,
            "confined_pitch_standard_error": error(completed * (end[:, 4] - mean_pitch))
            / survival
            if survival and args.copies > 1
            else None,
        }

    report["population_moments"] = {
        "SIMPLE": moments(simple, simple_codes == 0, simple_codes != 0),
        "ASCOT5": moments(ascot, ascot_finished, ascot_lost)
        if not ascot_unresolved.any()
        else "Not computed: unresolved or aborted ASCOT markers",
    }
    np.savetxt(
        output / "states.csv",
        np.column_stack((positions, simple, ascot)),
        delimiter=",",
        header="initial(R,phi,Z,E,pitch),SIMPLE(R,phi,Z,E,pitch),ASCOT(R,phi,Z,E,pitch)",
    )
    (output / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    with imas.DBEntry(str(output / "output.nc"), "w", dd_version="4.1.1") as entry:
        entry.put(result)
    print(json.dumps(report, indent=2))
    if ascot_unresolved.any():
        raise RuntimeError(
            f"{int(np.sum(ascot_unresolved))} ASCOT markers ended unresolved or "
            f"aborted ({ascot_states}); see ascot_native_errors in the report"
        )


if __name__ == "__main__":
    main()
