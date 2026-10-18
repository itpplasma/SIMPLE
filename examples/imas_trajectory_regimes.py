"""Regime-aware SIMPLE/ASCOT5 comparison of collisionless guiding-centre orbits.

Pointwise trajectory agreement is meaningful only while each code's own orbit
is insensitive to perturbations far below the comparison tolerance. Each code
is therefore rerun with a displaced start (``perturbation`` metres in R) and
with a different numerical resolution (SIMPLE relerr, ASCOT step halving).
The per-marker predictability horizon ``t_pred`` is the first sample at which
either self-divergence exceeds ``horizon_tolerance`` metres. Codes are compared
pointwise only for regular confined orbits and only for ``t <= t_pred``; lost,
near-separatrix and sensitive orbits enter statistics only.

ASCOT5 integrates relativistic guiding-centre equations, SIMPLE the
nonrelativistic ones. In a static field without E, relativistic guiding-centre
motion equals nonrelativistic motion with mass gamma*m at the same speed, so a
"relativistic-equivalent" SIMPLE run (n_d = A*gamma, speed set to ASCOT's)
isolates field/integration differences from the model difference.

Statistics use an independent random marker cohort (seeded) with paired tests:
loss fraction vs time (Wilson 95 % intervals, exact McNemar), loss-time and
strike-point distributions (two-sample Kolmogorov-Smirnov, Anderson-Darling)
and confined end-state moments (paired standard errors).
"""

import json
from pathlib import Path

import numpy as np

EV = 1.602176634e-19  # J, used by a5py/unyt for eV inputs
LIGHT = 299792458.0
PROTON = 1.6726e-27  # Benchmark marker mass unit in both codes
SIMPLE_EV_CGS = 1.6022e-12  # SIMPLE src/util.F90 ``ev``
SIMPLE_PROTON_CGS = 1.6726e-24  # SIMPLE src/util.F90 ``p_mass``
SIMPLE_EXIT = {0: "confined", 1: "LCFS", 2: "wall"}


def ascot_speed(energy, mass_number):
    """Relativistic speed [m/s] and gamma of a marker of kinetic energy [eV]."""
    gamma = 1 + energy * EV / (mass_number * PROTON * LIGHT**2)
    return LIGHT * np.sqrt(1 - 1 / gamma**2), gamma


def simple_speed(energy, mass_number):
    """SIMPLE's nonrelativistic speed [m/s] with its CGS constants."""
    return 1e-2 * np.sqrt(
        2 * energy * SIMPLE_EV_CGS / (mass_number * SIMPLE_PROTON_CGS)
    )


class CommonField:
    """IDS psi/B_phi cubic splines (ASCOT's 2D field interpolant) for diagnostics.

    DD4 (COCOS 17): B_R = -psi_Z/(2 pi R), B_Z = psi_R/(2 pi R) and the
    covariant vector potential is psi/(2 pi), so P_phi = m R v_par b_phi
    + q psi/(2 pi).
    """

    def __init__(self, equilibrium):
        from scipy.interpolate import RectBivariateSpline

        t = equilibrium.time_slice[0]
        p = t.profiles_2d[0]
        self.psi = RectBivariateSpline(p.grid.dim1, p.grid.dim2, p.psi)
        self.bphi = RectBivariateSpline(p.grid.dim1, p.grid.dim2, p.b_field_phi)
        self.axis = (
            float(t.global_quantities.magnetic_axis.r),
            float(t.global_quantities.magnetic_axis.z),
        )
        self.psi_axis = float(t.global_quantities.psi_axis)
        self.psi_edge = float(t.global_quantities.psi_boundary)

    def evaluate(self, r, z):
        psi = self.psi.ev(r, z)
        br = -self.psi.ev(r, z, dy=1) / (2 * np.pi * r)
        bz = self.psi.ev(r, z, dx=1) / (2 * np.pi * r)
        bphi = self.bphi.ev(r, z)
        return psi, bphi, np.sqrt(br**2 + bz**2 + bphi**2)

    def rho_pol(self, r, z):
        norm = (self.psi.ev(r, z) - self.psi_axis) / (self.psi_edge - self.psi_axis)
        return np.sqrt(np.clip(norm, 0, None))


def _simple_settings(mass_number, charge, duration, samples, energy, speed, relerr):
    import pysimple

    settings = dict(
        contr_pp=-1e30,
        notrace_passing=0,
        n_d=mass_number,
        n_e=charge,
        trace_time=duration,
        ntestpart=1,
        facE_al=1.0,
        orbit_model=0,
        integmode=pysimple.RK45,
        ntimstep=samples,
        deterministic=True,
        ran_seed=42,
    )
    if speed is not None:
        # The actor's state is v/v0 = sqrt(E/3.5 MeV) with
        # v0**2 = 2 (3.5 MeV/facE_al) ev/(n_d m_p) (src/params.f90).
        v0 = 100 * speed * np.sqrt(3.5e6 / energy)
        settings["facE_al"] = (
            2 * 3.5e6 * SIMPLE_EV_CGS / (mass_number * SIMPLE_PROTON_CGS * v0**2)
        )
    if relerr is not None:
        settings["relerr"] = relerr
    return settings


def simple_orbits(
    field,
    positions,
    mass_number,
    charge,
    duration,
    samples,
    *,
    speed=None,
    relerr=None,
):
    """Sampled SIMPLE RK45 orbits, one native ``trace_orbit`` per marker.

    ``mass_number`` may be non-integer (relativistic-equivalent mass).
    ``speed`` [m/s] overrides SIMPLE's energy-to-speed conversion.
    """
    import pysimple

    energy = float(positions[0, 3])
    if not np.allclose(positions[:, 3], energy, rtol=0, atol=0):
        raise ValueError("Regime analysis expects one marker energy")
    pysimple.init(
        field,
        **_simple_settings(
            mass_number, charge, duration, samples, energy, speed, relerr
        ),
    )
    v0 = float(pysimple.params.v0) * 1e-2
    states = np.empty((5, len(positions)), order="F")
    states[:3] = pysimple.cylindrical_to_reference(positions[:, :3].T)
    states[3], states[4] = np.sqrt(positions[:, 3] / 3.5e6), positions[:, 4]
    orbits = []
    for state in states.T:
        result = pysimple.trace_orbit(state, return_trajectory=True)
        code = int(np.asarray(pysimple.params.orbit_exit_code).ravel()[0])
        trajectory, times = result["trajectory"], result["times"]
        valid = np.all(np.isfinite(trajectory), axis=0) & np.isfinite(times)
        n = len(valid) if valid.all() else int(np.argmin(valid))
        # A loss sample sits at the exit time, off the shared output grid.
        grid = times[:n] * (samples - 1) / duration
        n = int(np.sum(np.abs(grid - np.rint(grid)) < 1e-6))
        r, phi, z = (
            pysimple.reference_to_cylindrical(trajectory[:3, :n])
            if n
            else np.zeros((3, 0))
        )
        orbit = {
            "t": times[:n],
            "r": r,
            "phi": phi,
            "z": z,
            "pitch": trajectory[4, :n],
            "speed": v0 * trajectory[3, :n],
            "status": SIMPLE_EXIT.get(code, f"numerical({code})"),
            "loss_time": result["loss_time"] if code in (1, 2) else None,
            "strike": None,
        }
        if code in (1, 2):
            end = np.asarray(result["final_position"])
            orbit["strike"] = pysimple.reference_to_cylindrical(end[:3]).tolist()
        orbits.append(orbit)
    return orbits, v0


def simple_cohort(field, positions, mass_number, charge, duration, *, speed=None):
    """SIMPLE end states (as the IDS actor's ``trace_parallel``) of a cohort."""
    import pysimple

    energy = float(positions[0, 3])
    settings = _simple_settings(
        mass_number, charge, duration, 1000, energy, speed, None
    )
    settings["ntestpart"] = len(positions)
    pysimple.init(field, **settings)
    states = np.empty((5, len(positions)), order="F")
    states[:3] = pysimple.cylindrical_to_reference(positions[:, :3].T)
    states[3], states[4] = np.sqrt(positions[:, 3] / 3.5e6), positions[:, 4]
    result = pysimple.trace_parallel(states)
    codes = np.asarray(pysimple.params.orbit_exit_code).copy()
    end = result["final_positions"]
    cylindrical = np.full((3, len(positions)), np.nan)
    ok = np.all(np.isfinite(end), axis=0) & np.isin(codes, [0, 1, 2])
    if ok.any():
        cylindrical[:, ok] = pysimple.reference_to_cylindrical(end[:3, ok])
    return {
        "status": [SIMPLE_EXIT.get(int(c), f"numerical({int(c)})") for c in codes],
        "lost": np.isin(codes, [1, 2]),
        "resolved": np.isin(codes, [0, 1, 2]),
        "time": np.where(np.isin(codes, [1, 2]), result["loss_times"], duration),
        "r": cylindrical[0],
        "z": cylindrical[2],
        "pitch": end[4].copy(),
    }


def ascot_orbits(
    equilibrium,
    positions,
    mass_number,
    charge,
    path,
    duration,
    step,
    samples=None,
):
    """ASCOT5 fixed-step relativistic GC orbits sampled on SIMPLE's time grid.

    Inputs mirror ``ascot_run`` in imas_benchmark.py without collisions. The
    recorded orbit is kept only at multiples of duration/(samples-1); the
    end-state supplies loss time and strike point.
    """
    from a5py import Ascot
    from a5py.ascot5io.options import Opt
    from a5py.ascot5io.state import State
    from pysimple.imas import Markers, write_markers
    import unyt

    markers = write_markers(
        [
            Markers(
                positions,
                np.ones(len(positions)),
                mass_number,
                charge,
                charge,
                rate=True,
            )
        ]
    )
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
        a.data.create_input(name)
    a.simulation_initinputs()
    options = Opt.get_default()
    options.update(
        SIM_MODE=2,
        ENABLE_ORBIT_FOLLOWING=1,
        ENABLE_ADAPTIVE=0,
        FIXEDSTEP_USE_USERDEFINED=1,
        FIXEDSTEP_USERDEFINED=step,
        ENDCOND_SIMTIMELIM=1,
        ENDCOND_LIM_SIMTIME=duration - step / 2,
        ENDCOND_RHOLIM=1,
        ENDCOND_MAX_RHO=1.0,
        ENDCOND_MIN_RHO=0.0,
    )
    if samples is not None:
        interval = duration / (samples - 1)
        every = round(interval / step)
        if every < 1 or not np.isclose(every * step, interval, rtol=1e-9, atol=0):
            raise ValueError("Orbit sampling interval must be a multiple of the step")
        # Recording happens once the elapsed time reaches the interval; half a
        # step margin keeps floating-point mileage from skipping a sample.
        options.update(
            ENABLE_ORBITWRITE=1,
            ORBITWRITE_MODE=1,
            ORBITWRITE_NPOINT=samples + 2,
            ORBITWRITE_INTERVAL=(every - 0.5) * step,
        )
    a.simulation_initoptions(**options)
    m["mass"] = np.full(m["n"], mass_number * PROTON) * unyt.kg
    m["zeta"] = np.zeros(m["n"]) * unyt.rad
    a.simulation_initmarkers(**m)
    run = a.simulation_run(printsummary=False)
    r, phi, z, pitch, condition, end_time, errors = run.getstate(
        "r", "phi", "z", "pitch", "endcond", "time", "errormsg", state="end", mode="gc"
    )
    codes = np.asarray(condition).astype(int)
    lost = np.array([State.endcond_check(int(c), "RHOMAX") for c in codes])
    finished = np.array([State.endcond_check(int(c), "TLIM") for c in codes])
    resolved = (lost != finished) & (np.asarray(errors).astype(int) == 0)
    end = {
        "status": [
            State.endcond_tostring(int(c)) or f"UNKNOWN({int(c)})" for c in codes
        ],
        "lost": lost & resolved,
        "resolved": resolved,
        "time": end_time.to_value("s").copy(),
        "r": r.to_value("m").copy(),
        "phi": phi.to_value("rad").copy(),
        "z": z.to_value("m").copy(),
        "pitch": pitch.to_value("1").copy(),
    }
    orbits = None
    if samples is not None:
        ids, t, ro, po, zo, pi, ek = run.getorbit(
            "ids", "mileage", "r", "phi", "z", "pitch", "ekin"
        )
        ids, t = np.asarray(ids).astype(int), t.to_value("s")
        index = t / interval
        on_grid = np.abs(index - np.rint(index)) < 1e-6
        orbits = []
        for i in range(len(positions)):
            select = np.flatnonzero((ids == i + 1) & on_grid)
            select = select[np.argsort(t[select])]
            orbits.append(
                {
                    "t": t[select],
                    "r": ro.to_value("m")[select],
                    "phi": po.to_value("rad")[select],
                    "z": zo.to_value("m")[select],
                    "pitch": pi.to_value("1")[select],
                    "energy": ek.to_value("eV")[select],
                    "status": end["status"][i] if resolved[i] else "unresolved",
                    "loss_time": float(end["time"][i]) if end["lost"][i] else None,
                    "strike": [end["r"][i], end["phi"][i], end["z"][i]]
                    if end["lost"][i]
                    else None,
                }
            )
    a.simulation_free(inputs=True, markers=True, diagnostics=True)
    return orbits, end


def _cartesian(orbit, n):
    r, phi = orbit["r"][:n], orbit["phi"][:n]
    return np.column_stack((r * np.cos(phi), r * np.sin(phi), orbit["z"][:n]))


def _common(a, b):
    """Leading samples present in both orbits on the shared time grid."""
    n = min(len(a["t"]), len(b["t"]))
    if n and not np.allclose(a["t"][:n], b["t"][:n], rtol=0, atol=1e-6 * a["t"][-1]):
        raise ValueError("Orbit samples are not on a common time grid")
    return n


def distance(a, b):
    n = _common(a, b)
    return a["t"][:n], np.linalg.norm(_cartesian(a, n) - _cartesian(b, n), axis=1)


def _crossings(t, values):
    """Interpolated zero crossings: (time, direction, fraction, index)."""
    sign = np.sign(values)
    k = np.flatnonzero(sign[:-1] * sign[1:] < 0)
    fraction = values[k] / (values[k] - values[k + 1])
    return t[k] + fraction * (t[k + 1] - t[k]), sign[k + 1], fraction, k


def orbit_averages(orbit, field):
    """Bounce/transit period, toroidal precession, turning points and width."""
    t, pitch = orbit["t"], orbit["pitch"]
    if len(t) < 8:
        return {"class": "unresolved_short"}
    phi = np.unwrap(orbit["phi"])
    rho = field.rho_pol(orbit["r"], orbit["z"])
    times, direction, fraction, k = _crossings(t, pitch)
    result = {"radial_width_rho_pol": float(np.ptp(rho))}
    if len(times) == 0:
        theta = np.unwrap(
            np.arctan2(orbit["z"] - field.axis[1], orbit["r"] - field.axis[0])
        )
        result["class"] = "passing"
        result["min_abs_pitch"] = float(np.min(np.abs(pitch)))
        if abs(theta[-1] - theta[0]) >= 4 * np.pi:
            result["transit_period_s"] = float(
                2 * np.pi / abs(np.polyfit(t, theta, 1)[0])
            )
            result["toroidal_frequency_rad_s"] = float(np.polyfit(t, phi, 1)[0])
        return result
    result["class"] = "trapped"
    same = times[direction > 0]
    if len(same) >= 2:
        phase = phi[k] + fraction * (phi[k + 1] - phi[k])
        phase = phase[direction > 0]
        result["bounce_period_s"] = float(np.mean(np.diff(same)))
        result["bounce_period_spread_s"] = float(np.ptp(np.diff(same)))
        result["precession_rad_per_bounce"] = float(np.mean(np.diff(phase)))
    for name, sign in (("turning_point_up", 1), ("turning_point_down", -1)):
        select = direction == sign
        if select.any():
            rk = orbit["r"][k] + fraction * (orbit["r"][k + 1] - orbit["r"][k])
            zk = orbit["z"][k] + fraction * (orbit["z"][k + 1] - orbit["z"][k])
            result[name] = [float(np.mean(rk[select])), float(np.mean(zk[select]))]
    return result


def invariant_drift(orbit, field, mass, charge, speed, gamma=1.0):
    """Maximum drift of mu and P_phi along an orbit, from the common field.

    Momentum is gamma*m*v (gamma = 1 for SIMPLE). P_phi is normalized by
    q*|psi_edge - psi_axis|/(2 pi).
    """
    if len(orbit["t"]) < 2:
        return None
    psi, bphi, b = field.evaluate(orbit["r"], orbit["z"])
    v = orbit["speed"] if "speed" in orbit else speed
    p = gamma * mass * v
    mu = p**2 * (1 - orbit["pitch"] ** 2) / (2 * mass * b)
    ptor = orbit["r"] * p * orbit["pitch"] * bphi / b + charge * psi / (2 * np.pi)
    scale = charge * abs(field.psi_edge - field.psi_axis) / (2 * np.pi)
    out = {"P_phi_relative_drift": float(np.max(np.abs(ptor - ptor[0])) / scale)}
    if mu[0] > 1e-12 * mass * np.max(v) ** 2:
        out["mu_relative_drift"] = float(np.max(np.abs(mu / mu[0] - 1)))
    if "energy" in orbit:
        out["energy_relative_drift"] = float(
            np.max(np.abs(orbit["energy"] / orbit["energy"][0] - 1))
        )
    return out


def _horizon(t, d, tolerance):
    above = np.flatnonzero(d > tolerance)
    return (float(t[above[0]]), True) if len(above) else (float(t[-1]), False)


def _growth(t, d, delta):
    """Finite-time Lyapunov estimate and log-log growth exponent."""
    ok = (t > 0) & (d > 0)
    if ok.sum() < 4:
        return None, None
    t, d = t[ok], d[ok]
    lyapunov = float(np.log(d[-1] / delta) / t[-1])
    late = t >= t[-1] / 10
    exponent = (
        float(np.polyfit(np.log(t[late]), np.log(d[late]), 1)[0])
        if late.sum() >= 4
        else None
    )
    return lyapunov, exponent


def _relative(a, b):
    if a is None or b is None:
        return None
    return float(abs(a - b) / max(abs(b), 1e-300))


def compare_averages(a, b):
    """Relative differences of matching orbit-averaged quantities."""
    out = {}
    for key in [
        "bounce_period_s",
        "transit_period_s",
        "toroidal_frequency_rad_s",
        "precession_rad_per_bounce",
        "radial_width_rho_pol",
    ]:
        if key in a and key in b:
            out[key] = _relative(a[key], b[key])
    for key in ["turning_point_up", "turning_point_down"]:
        if key in a and key in b:
            out[key + "_distance_m"] = float(np.hypot(*np.subtract(a[key], b[key])))
    return out


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return None
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [float(p), float(centre - half), float(centre + half)]


def mcnemar(a, b):
    """Exact two-sided McNemar test for paired binary outcomes."""
    from scipy.stats import binom

    only_a, only_b = int(np.sum(a & ~b)), int(np.sum(~a & b))
    n = only_a + only_b
    p = 1.0 if n == 0 else min(1.0, 2 * binom.cdf(min(only_a, only_b), n, 0.5))
    return {"only_first": only_a, "only_second": only_b, "p_value": float(p)}


def two_sample(x, y):
    from scipy.stats import anderson_ksamp, ks_2samp

    if len(x) < 5 or len(y) < 5:
        return {"n": [len(x), len(y)], "note": "fewer than 5 samples"}
    ks = ks_2samp(x, y)
    out = {
        "n": [len(x), len(y)],
        "ks_statistic": float(ks.statistic),
        "ks_p_value": float(ks.pvalue),
    }
    try:
        ad = anderson_ksamp([x, y])
        # SciPy caps the interpolated significance to [0.001, 0.25].
        out["anderson_darling_p_value_capped"] = float(ad.pvalue)
    except ValueError as error:
        out["anderson_darling"] = str(error)
    return out


def paired_moment(x, y):
    d = x - y
    if len(d) < 2:
        return None
    se = float(np.std(d, ddof=1) / np.sqrt(len(d)))
    return {
        "first": float(np.mean(x)),
        "second": float(np.mean(y)),
        "difference": float(np.mean(d)),
        "standard_error": se,
        "within_3_se": bool(abs(np.mean(d)) <= 3 * se)
        if se > 0
        else bool(np.mean(d) == 0),
    }


def statistics(first, second, duration, field, labels):
    """Paired cohort statistics; unresolved markers are excluded, not lost."""
    ok = first["resolved"] & second["resolved"]
    n = int(ok.sum())
    a_lost, b_lost = first["lost"][ok], second["lost"][ok]
    checkpoints = duration * np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
    curve = []
    for t in checkpoints:
        ka = int(np.sum(a_lost & (first["time"][ok] <= t)))
        kb = int(np.sum(b_lost & (second["time"][ok] <= t)))
        curve.append(
            {"time_s": float(t), labels[0]: wilson(ka, n), labels[1]: wilson(kb, n)}
        )

    def strike(end, lost):
        return np.arctan2(
            end["z"][ok][lost] - field.axis[1], end["r"][ok][lost] - field.axis[0]
        )

    both = ~a_lost & ~b_lost
    moments = {}
    for name, fa, fb in [
        ("pitch", first["pitch"][ok][both], second["pitch"][ok][both]),
        (
            "pitch_squared",
            first["pitch"][ok][both] ** 2,
            second["pitch"][ok][both] ** 2,
        ),
        (
            "rho_pol",
            field.rho_pol(first["r"][ok][both], first["z"][ok][both]),
            field.rho_pol(second["r"][ok][both], second["z"][ok][both]),
        ),
    ]:
        moments[name] = paired_moment(fa, fb)
    tests = {
        "markers_compared": n,
        "unresolved_excluded": {
            labels[0]: int(np.sum(~first["resolved"])),
            labels[1]: int(np.sum(~second["resolved"])),
        },
        "loss_fraction_vs_time_wilson95": curve,
        "lost_by_end_mcnemar": mcnemar(a_lost, b_lost),
        "loss_time_distribution": two_sample(
            first["time"][ok][a_lost], second["time"][ok][b_lost]
        ),
        "strike_poloidal_angle_distribution": two_sample(
            strike(first, a_lost), strike(second, b_lost)
        ),
        "confined_in_both_paired_moments": moments,
    }
    passed = [tests["lost_by_end_mcnemar"]["p_value"] >= 0.01]
    for key in ["loss_time_distribution", "strike_poloidal_angle_distribution"]:
        if "ks_p_value" in tests[key]:
            passed.append(tests[key]["ks_p_value"] >= 0.01)
    passed += [m["within_3_se"] for m in moments.values() if m is not None]
    tests["accepted"] = bool(all(passed))
    tests["acceptance"] = "McNemar and KS p >= 0.01; paired moments within 3 SE"
    return tests


def _class_signature(averages):
    return averages.get("class")


def marker_regime(record, duration):
    """Classify one marker from both codes' nominal and perturbed orbits."""
    if record["unresolved"]:
        return "unresolved"
    if record["lost"]:
        return "lost"
    classes = set(record["classes"].values())
    if len(classes) > 1:
        return "near_separatrix"
    if classes == {"passing"} and record["min_abs_pitch"] < 0.05:
        return "near_separatrix"
    if record["t_pred_s"] < duration * (1 - 1e-9):
        exponent = record["growth_exponent"]
        return "stochastic" if exponent is not None and exponent > 2 else "sensitive"
    return f"regular_{classes.pop()}"


def compare(
    equilibrium,
    field,
    positions,
    mass_number,
    charge,
    *,
    duration,
    step,
    output,
    samples=1001,
    perturbation=1e-6,
    horizon_tolerance=1e-3,
    relative_tolerance=1e-4,
    absolute_tolerance=1e-4,
    stat_markers=0,
    seed=12345,
):
    """Run all variants, classify markers and return the regime report."""
    output = Path(output)
    common = CommonField(equilibrium)
    energy = float(positions[0, 3])
    mass, q = mass_number * PROTON, charge * EV
    v_ascot, gamma = ascot_speed(energy, mass_number)
    perturbed = positions.copy()
    perturbed[:, 0] += perturbation

    simple = {}
    simple["nominal"], v_simple = simple_orbits(
        field, positions, mass_number, charge, duration, samples
    )
    simple["perturbed"], _ = simple_orbits(
        field, perturbed, mass_number, charge, duration, samples
    )
    simple["relerr_1e-10"], _ = simple_orbits(
        field, positions, mass_number, charge, duration, samples, relerr=1e-10
    )
    simple["relativistic_equivalent"], v_equivalent = simple_orbits(
        field, positions, mass_number * gamma, charge, duration, samples, speed=v_ascot
    )
    both, _ = ascot_orbits(
        equilibrium,
        np.vstack((positions, perturbed)),
        mass_number,
        charge,
        output / "ascot-orbits.h5",
        duration,
        step,
        samples,
    )
    half, _ = ascot_orbits(
        equilibrium,
        positions,
        mass_number,
        charge,
        output / "ascot-orbits-half-step.h5",
        duration,
        step / 2,
        samples,
    )
    ascot = {
        "nominal": both[: len(positions)],
        "perturbed": both[len(positions) :],
        "half_step": half,
    }

    markers = []
    for i, start in enumerate(positions):
        s = {k: v[i] for k, v in simple.items()}
        a = {k: v[i] for k, v in ascot.items()}
        unresolved = any(
            o["status"].startswith(("numerical", "unresolved"))
            for o in [*s.values(), *a.values()]
        )
        lost = any(o["loss_time"] is not None for o in [s["nominal"], a["nominal"]])
        horizons, growth = {}, {}
        for code, nominal, ic, numeric in [
            ("SIMPLE", s["nominal"], s["perturbed"], s["relerr_1e-10"]),
            ("ASCOT5", a["nominal"], a["perturbed"], a["half_step"]),
        ]:
            t, d_ic = distance(nominal, ic)
            _, d_num = distance(nominal, numeric)
            if len(t) == 0:
                horizons[code] = (0.0, True)
                continue
            n = min(len(d_ic), len(d_num))
            horizons[code] = _horizon(
                t[:n], np.maximum(d_ic[:n], d_num[:n]), horizon_tolerance
            )
            growth[code] = {
                "max_ic_divergence_m": float(np.max(d_ic)),
                "max_numerical_divergence_m": float(np.max(d_num)),
                "amplification": float(np.max(d_ic) / perturbation),
            }
            (
                growth[code]["finite_time_lyapunov_per_s"],
                growth[code]["growth_exponent"],
            ) = _growth(t, d_ic, perturbation)
        t_pred = min(h[0] for h in horizons.values())
        averages = {
            "SIMPLE": orbit_averages(s["nominal"], common),
            "SIMPLE_perturbed": orbit_averages(s["perturbed"], common),
            "SIMPLE_relativistic_equivalent": orbit_averages(
                s["relativistic_equivalent"], common
            ),
            "ASCOT5": orbit_averages(a["nominal"], common),
            "ASCOT5_perturbed": orbit_averages(a["perturbed"], common),
        }
        exponents = [
            g["growth_exponent"]
            for g in growth.values()
            if g.get("growth_exponent") is not None
        ]
        record = {
            "marker": i,
            "start": start.tolist(),
            "unresolved": unresolved,
            "lost": lost,
            "status": {
                "SIMPLE": s["nominal"]["status"],
                "ASCOT5": a["nominal"]["status"],
            },
            "loss_time_s": {
                "SIMPLE": s["nominal"]["loss_time"],
                "ASCOT5": a["nominal"]["loss_time"],
            },
            "classes": {k: _class_signature(v) for k, v in averages.items()},
            "min_abs_pitch": min(
                v.get("min_abs_pitch", 1.0) for v in averages.values()
            ),
            "t_pred_s": t_pred,
            "horizon_reached": {k: h[1] for k, h in horizons.items()},
            "self_divergence": growth,
            "growth_exponent": max(exponents) if exponents else None,
        }
        record["regime"] = marker_regime(record, duration)
        # Pointwise comparison only inside the horizon and before any loss.
        window = {}
        for name, simple_orbit in [
            ("nominal", s["nominal"]),
            ("relativistic_equivalent", s["relativistic_equivalent"]),
        ]:
            t, d = distance(simple_orbit, a["nominal"])
            inside = t <= t_pred
            if inside.any():
                path = float(np.mean(simple_orbit["speed"])) * float(t[inside][-1])
                window[name] = {
                    "window_s": float(t[inside][-1]),
                    "max_distance_m": float(np.max(d[inside])),
                    "relative_to_path": float(np.max(d[inside]) / max(path, 1e-300)),
                }
        if "nominal" in window:
            sn = s["nominal"]
            n = int(np.sum(sn["t"] <= window["nominal"]["window_s"]))
            # Leading-order prediction: parallel phase offset from the speed
            # difference, integrated along the SIMPLE orbit.
            window["nominal"]["speed_offset_prediction_m"] = (
                float(
                    abs(v_simple - v_ascot)
                    * np.trapezoid(np.abs(sn["pitch"][:n]), sn["t"][:n])
                )
                if n > 1
                else None
            )
        record["pointwise"] = window
        start_gap = (
            float(
                np.linalg.norm(
                    _cartesian(s["nominal"], 1) - _cartesian(a["nominal"], 1)
                )
            )
            if len(s["nominal"]["t"]) and len(a["nominal"]["t"])
            else None
        )
        record["start_position_gap_m"] = start_gap
        if record["regime"].startswith("regular"):
            eq_window = window.get("relativistic_equivalent")
            record["deterministic_accepted"] = bool(
                eq_window is not None
                and eq_window["max_distance_m"]
                <= max(
                    absolute_tolerance,
                    relative_tolerance * v_ascot * eq_window["window_s"],
                )
            )
            record["orbit_average_relative_differences"] = {
                "SIMPLE_vs_ASCOT5": compare_averages(
                    averages["SIMPLE"], averages["ASCOT5"]
                ),
                "relativistic_equivalent_vs_ASCOT5": compare_averages(
                    averages["SIMPLE_relativistic_equivalent"], averages["ASCOT5"]
                ),
                "ASCOT5_perturbed_vs_ASCOT5": compare_averages(
                    averages["ASCOT5_perturbed"], averages["ASCOT5"]
                ),
            }
            record["orbit_averages"] = {
                "SIMPLE": averages["SIMPLE"],
                "ASCOT5": averages["ASCOT5"],
                "SIMPLE_relativistic_equivalent": averages[
                    "SIMPLE_relativistic_equivalent"
                ],
            }
            record["invariant_drift"] = {
                "SIMPLE": invariant_drift(s["nominal"], common, mass, q, v_simple),
                "SIMPLE_relativistic_equivalent": invariant_drift(
                    s["relativistic_equivalent"], common, mass * gamma, q, v_equivalent
                ),
                "ASCOT5": invariant_drift(
                    a["nominal"], common, mass, q, v_ascot, gamma
                ),
            }
        markers.append(record)

    regimes = {}
    for record in markers:
        regimes[record["regime"]] = regimes.get(record["regime"], 0) + 1
    regular = [m for m in markers if m["regime"].startswith("regular")]

    def worst(path):
        values = []
        for m in regular:
            value = m
            for key in path:
                value = value.get(key) if isinstance(value, dict) else None
            if value is not None:
                values.append(value)
        return float(np.max(values)) if values else None

    summary = {
        "regime_counts": regimes,
        "t_pred_s_quantiles_all": np.quantile(
            [m["t_pred_s"] for m in markers], [0, 0.25, 0.5, 0.75, 1]
        ).tolist(),
        "regular_markers": len(regular),
        "regular_deterministic_accepted": int(
            sum(m["deterministic_accepted"] for m in regular)
        ),
        "regular_max_distance_m": {
            "nominal": worst(["pointwise", "nominal", "max_distance_m"]),
            "relativistic_equivalent": worst(
                ["pointwise", "relativistic_equivalent", "max_distance_m"]
            ),
        },
        "regular_max_relative_to_path_relativistic_equivalent": worst(
            ["pointwise", "relativistic_equivalent", "relative_to_path"]
        ),
        "regular_max_orbit_average_relative_difference": {
            key: {
                pair: worst(["orbit_average_relative_differences", pair, key])
                for pair in [
                    "SIMPLE_vs_ASCOT5",
                    "relativistic_equivalent_vs_ASCOT5",
                    "ASCOT5_perturbed_vs_ASCOT5",
                ]
            }
            for key in [
                "bounce_period_s",
                "transit_period_s",
                "toroidal_frequency_rad_s",
                "precession_rad_per_bounce",
                "radial_width_rho_pol",
                "turning_point_up_distance_m",
                "turning_point_down_distance_m",
            ]
        },
        "regular_max_invariant_drift": {
            code: {
                key: worst(["invariant_drift", code, key])
                for key in [
                    "P_phi_relative_drift",
                    "mu_relative_drift",
                    "energy_relative_drift",
                ]
            }
            for code in ["SIMPLE", "SIMPLE_relativistic_equivalent", "ASCOT5"]
        },
        "max_start_position_gap_m": float(
            np.max(
                [
                    m["start_position_gap_m"]
                    for m in markers
                    if m["start_position_gap_m"] is not None
                ]
            )
        ),
    }
    report = {
        "method": __doc__.strip(),
        "samples": samples,
        "duration_s": duration,
        "ascot_step_s": step,
        "perturbation_m": perturbation,
        "horizon_tolerance_m": horizon_tolerance,
        "deterministic_tolerance": {
            "relative_to_path": relative_tolerance,
            "absolute_m": absolute_tolerance,
        },
        "speeds_m_s": {
            "SIMPLE": v_simple,
            "ASCOT5": v_ascot,
            "SIMPLE_relativistic_equivalent": v_equivalent,
        },
        "ascot_gamma": float(gamma),
        "relative_speed_offset_SIMPLE_minus_ASCOT5": float(v_simple / v_ascot - 1),
        "summary": summary,
        "markers": markers,
    }
    if stat_markers:
        rng = np.random.default_rng(seed)
        t = equilibrium.time_slice[0]
        minor = max(t.boundary.outline.r) - common.axis[0]
        cohort = np.column_stack(
            (
                common.axis[0] + minor * rng.uniform(0.05, 0.97, stat_markers),
                np.full(stat_markers, 0.3),
                np.full(stat_markers, common.axis[1]),
                np.full(stat_markers, energy),
                rng.uniform(-1, 1, stat_markers),
            )
        )
        s_nominal = simple_cohort(field, cohort, mass_number, charge, duration)
        s_equivalent = simple_cohort(
            field, cohort, mass_number * gamma, charge, duration, speed=v_ascot
        )
        _, a_end = ascot_orbits(
            equilibrium,
            cohort,
            mass_number,
            charge,
            output / "ascot-cohort.h5",
            duration,
            step,
        )
        report["statistics"] = {
            "cohort": {
                "markers": stat_markers,
                "seed": seed,
                "sampling": "uniform outboard-midplane minor radius 0.05-0.97 "
                "and pitch in [-1, 1]",
            },
            "SIMPLE_vs_ASCOT5": statistics(
                s_nominal, a_end, duration, common, ["SIMPLE", "ASCOT5"]
            ),
            "relativistic_equivalent_vs_ASCOT5": statistics(
                s_equivalent,
                a_end,
                duration,
                common,
                ["SIMPLE_relativistic_equivalent", "ASCOT5"],
            ),
        }
    return report


def write(report, output):
    path = Path(output) / "regimes.json"
    path.write_text(json.dumps(report, indent=1, allow_nan=False) + "\n")
    return path
