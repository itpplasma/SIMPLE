"""SIMPLE's frozen-equilibrium IMAS cohort actor and database entry point."""

from pathlib import Path
import tempfile
import numpy as np

from ._common import factory, ids_version, initialize, require_imas
from .equilibrium import equilibrium_to_chartmap
from .markers import fill_marker_slice, fill_species, read_markers
from .plasma import plasma_profiles
from .wall import wall_to_stl

EV_TO_J = 1.602176634e-19


def run(
    markers,
    *,
    equilibrium=None,
    field=None,
    core_profiles=None,
    wall=None,
    duration,
    source_interval=None,
    time=None,
    workdir=None,
    conversion=None,
    **settings,
):
    """Follow guiding-centre ion cohorts and return a distributions IDS.

    Supply either an axisymmetric equilibrium IDS or an existing SIMPLE field
    file (VMEC/chartmap). Each source weight [particles/s] is multiplied by the
    explicit ``source_interval`` [s]. This is a finite cohort calculation, not
    a continuous-source steady-state distribution or a residence histogram.

    Output contains initial/final particle snapshots and global particle/energy
    totals. Lost markers have zero final weight. Private coordinate identifiers
    -1/-2 carry the event time [s] and the native exit condition; numerical or
    skipped markers raise RuntimeError, never silently count as physical losses.
    The equilibrium and profiles are frozen at the selected initial time.
    """
    import pysimple

    if not np.isfinite(duration) or duration <= 0:
        raise ValueError("duration must be positive and finite")
    if (equilibrium is None) == (field is None):
        raise ValueError("Supply exactly one equilibrium IDS or SIMPLE field file")
    populations = read_markers(markers, time=time)
    start_time = populations[0].time
    if populations[0].rate:
        if (
            source_interval is None
            or not np.isfinite(source_interval)
            or source_interval <= 0
        ):
            raise ValueError(
                "Source weights [particles/s] require positive source_interval [s]"
            )
    elif source_interval is not None:
        raise ValueError("source_interval applies only to distribution_sources weights")
    protected = {
        "n_d",
        "n_e",
        "trace_time",
        "ntestpart",
        "facE_al",
        "field_input",
        "coord_input",
        "netcdffile",
        "wall_input",
        "wall_units",
        "plasma_profiles",
        "orbit_model",
    }
    if protected.intersection(settings):
        raise ValueError(
            f"Actor owns settings {sorted(protected.intersection(settings))}"
        )
    if settings.get("swcoll", False) and core_profiles is None:
        raise ValueError("Collisions require a core_profiles IDS")
    if core_profiles is not None and field is not None:
        # An arbitrary preconverted field's toroidal-flux normalization must
        # be known before attaching independent IDS profiles.
        raise ValueError("IDS collisions require the matching equilibrium IDS")
    if core_profiles is not None and any(
        p.mass_number != 4 or p.charge != 2 or p.nuclear_charge != 2
        for p in populations
    ):
        raise ValueError(
            "SIMPLE's collision operator currently supports alpha markers (A=4, Z=2) only"
        )
    versions = {
        ids_version(i)
        for i in [markers, equilibrium, core_profiles, wall]
        if i is not None
    }
    if len(versions) != 1:
        raise ValueError("Convert input IDS to a common Data Dictionary version first")
    output = factory(ids_version(markers)).new("distributions")
    initialize(output, [start_time, start_time + duration])
    output.ids_properties.comment = (
        f"Frozen-equilibrium SIMPLE guiding-centre cohort; duration={duration} s; "
        f"source integration interval={source_interval} s; lost final weights are zero. "
        "This output is not a steady-state residence histogram."
    )
    output.distribution.resize(len(populations))
    with tempfile.TemporaryDirectory(prefix="simple-imas-", dir=workdir) as directory:
        directory = Path(directory)
        if equilibrium is not None:
            field = equilibrium_to_chartmap(
                equilibrium,
                directory / "equilibrium.nc",
                time=start_time,
                **(conversion or {}),
            )
        collision_table, collision_settings = None, {"swcoll": False}
        if core_profiles is not None:
            collision_table, collision_settings = plasma_profiles(
                core_profiles, time=start_time
            )
        wall_settings = {}
        if wall is not None:
            from pysimple._tracing import _api

            if not _api.simple_wall_supported():
                raise RuntimeError(
                    "IDS walls require SIMPLE built with SIMPLE_ENABLE_CGAL=ON"
                )
            path = wall_to_stl(wall, directory / "wall.stl", time=start_time)
            wall_settings = {"wall_input": str(path), "wall_units": "m"}
        for population, target in zip(populations, output.distribution):
            config = {
                "contr_pp": -1e30,
                "notrace_passing": 0,
                **collision_settings,
                **settings,
                **wall_settings,
                "n_d": population.mass_number,
                "n_e": population.charge,
                "trace_time": duration,
                "ntestpart": len(population.weights),
                "facE_al": 1.0,
                "orbit_model": 0,
            }
            pysimple.init(field, plasma_profiles=collision_table, **config)
            initial = population.positions
            states = np.empty((5, len(initial)), order="F")
            states[:3] = pysimple.cylindrical_to_reference(initial[:, :3].T)
            states[3] = np.sqrt(initial[:, 3] / 3.5e6)
            states[4] = initial[:, 4]
            result = pysimple.trace_parallel(states)
            codes = np.asarray(pysimple.params.orbit_exit_code).copy()
            if np.any(~np.isin(codes, [0, 1, 2])):
                failed = np.flatnonzero(~np.isin(codes, [0, 1, 2]))
                raise RuntimeError(
                    f"Untraced/numerically failed IDS markers: "
                    f"{list(zip(failed.tolist(), codes[failed].tolist()))}"
                )
            final = np.empty_like(initial)
            final[:, :3] = pysimple.reference_to_cylindrical(
                result["final_positions"][:3]
            ).T
            final[:, 3] = 3.5e6 * result["final_positions"][3] ** 2
            final[:, 4] = result["final_positions"][4]
            weights = population.weights * (source_interval if population.rate else 1)
            final_weights = weights * (codes == 0)
            target.gyro_type = 2
            fill_species(target.species, population)
            target.markers.resize(2)
            fill_marker_slice(target.markers[0], initial, weights, start_time)
            diagnostics = {
                "simple_event_time": (
                    -1,
                    "Absolute time [s] at termination",
                    start_time + result["loss_times"],
                ),
                "simple_exit_condition": (
                    -2,
                    "0 completed, 1 LCFS loss, 2 wall loss",
                    codes,
                ),
            }
            fill_marker_slice(
                target.markers[1],
                final,
                final_weights,
                start_time + duration,
                diagnostics=diagnostics,
            )
            target.global_quantities.resize(2)
            for q, t, position, weight in zip(
                target.global_quantities,
                output.time,
                [initial, final],
                [weights, final_weights],
            ):
                q.time = t
                q.particles_n = float(np.sum(weight))
                q.energy = float(np.dot(weight, position[:, 3]) * EV_TO_J)
    output.validate()
    return output


def run_entry(
    input_uri,
    output_uri,
    *,
    marker_ids="distribution_sources",
    occurrence=0,
    core_profiles=False,
    wall=False,
    dd_version="4.1.1",
    **options,
):
    """Read/write standard IMAS DBEntry URIs (or plain .nc paths).

    Output is created with mode ``x`` so an existing dataset is never replaced.
    ``occurrence`` selects the marker occurrence; other IDS use occurrence zero.
    """
    imas = require_imas()
    if marker_ids not in ("distribution_sources", "distributions"):
        raise ValueError("marker_ids must be distribution_sources or distributions")
    with imas.DBEntry(str(input_uri), "r", dd_version=dd_version) as entry:
        kwargs = {"equilibrium": entry.get("equilibrium")}
        if core_profiles:
            kwargs["core_profiles"] = entry.get("core_profiles")
        if wall:
            kwargs["wall"] = entry.get("wall")
        result = run(entry.get(marker_ids, occurrence), **kwargs, **options)
    with imas.DBEntry(str(output_uri), "x", dd_version=dd_version) as entry:
        entry.put(result)
    return result


def main(argv=None):
    """Command-line actor for an IMAS database or portable NetCDF entry."""
    import argparse
    import json

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("input_uri")
    parser.add_argument("output_uri")
    parser.add_argument("--duration", type=float, required=True)
    parser.add_argument("--source-interval", type=float)
    parser.add_argument("--time", type=float)
    parser.add_argument(
        "--marker-ids",
        default="distribution_sources",
        choices=["distribution_sources", "distributions"],
    )
    parser.add_argument("--dd-version", default="4.1.1")
    parser.add_argument("--core-profiles", action="store_true")
    parser.add_argument("--wall", action="store_true")
    parser.add_argument(
        "--settings", type=Path, help="JSON object of SIMPLE &config settings"
    )
    args = vars(parser.parse_args(argv))
    settings = args.pop("settings")
    if settings:
        args.update(json.loads(settings.read_text()))
    run_entry(**args)
