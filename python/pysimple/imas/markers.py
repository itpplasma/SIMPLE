"""Standard fast-ion marker mappings shared with ASCOT5's IMAS input."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from ._common import array, check_ids, factory, initialize, scalar, select_slice

# IMAS coordinate_identifier: metres, radians, metres, eV, v_parallel/v.
MARKER_COORDINATES = (4, 5, 3, 301, 403)
MARKER_NAMES = ("r", "phi", "z", "energy_kinetic", "pitch")


@dataclass
class Markers:
    """One ion population in physical cylindrical coordinates.

    ``positions`` has columns R [m], phi [rad], Z [m], kinetic energy [eV],
    pitch v_parallel/v. Weights are particles/s for distribution_sources,
    particles for distributions. Source weights need an explicit integration
    interval before they can describe a finite cohort in distributions.
    """

    positions: np.ndarray
    weights: np.ndarray
    mass_number: float
    nuclear_charge: int
    charge: int
    time: float = 0.0
    rate: bool = False

    def __post_init__(self):
        self.positions = array(self.positions, "marker positions", ndim=2)
        self.weights = array(self.weights, "marker weights", nonnegative=True)
        if self.positions.shape != (self.weights.size, 5):
            raise ValueError("marker positions must have shape (n, 5)")
        if np.any(self.positions[:, 0] <= 0) or np.any(self.positions[:, 3] <= 0):
            raise ValueError("marker R and kinetic energy must be positive")
        if np.any(np.abs(self.positions[:, 4]) > 1):
            raise ValueError("marker pitch must be in [-1, 1]")
        if not np.isfinite(self.mass_number) or self.mass_number <= 0:
            raise ValueError("mass_number must be positive")
        if not (1 <= self.charge <= self.nuclear_charge):
            raise ValueError("SIMPLE IDS markers require positively charged ions")
        if not np.isfinite(self.time):
            raise ValueError("marker time must be finite")


def _species(species):
    if int(species.type.index) not in (2, 3):
        raise ValueError("SIMPLE IDS input currently requires ion markers")
    ion = species.ion
    if len(ion.element) != 1:
        raise ValueError("Molecular/mixed-element marker species are unsupported")
    a = scalar(ion.element[0].a, "ion mass number")
    z = scalar(ion.element[0].z_n, "ion nuclear charge")
    if int(species.type.index) == 3:
        charge = scalar(ion.state.z_min, "ion charge state")
        if charge != scalar(ion.state.z_max, "ion charge state"):
            raise ValueError("Use one ion charge state per marker population")
    else:
        charge = scalar(ion.z_ion, "ion charge")
    if z != int(z) or charge != int(charge):
        raise ValueError("SIMPLE requires integer nuclear and ion charges")
    return a, int(z), int(charge)


def read_markers(ids, *, time=None):
    """Read all populations from distribution_sources or distributions.

    Uses coordinate *identifiers*, accepting reordered columns and additional
    coordinates. Actual-particle (gyro_type=1) data are rejected rather than
    silently interpreted as guiding centres. Time selection is exact.
    """
    name = ids.metadata.name
    if name not in ("distribution_sources", "distributions"):
        raise ValueError("Expected distribution_sources or distributions IDS")
    check_ids(ids, name)
    populations = ids.source if name == "distribution_sources" else ids.distribution
    result = []
    for population in populations:
        if not len(population.markers):
            raise ValueError("Each population must provide particle markers")
        if int(population.gyro_type) != 2:
            raise ValueError("SIMPLE's IDS actor requires gyro_type=2 guiding centres")
        markers = select_slice(population.markers, time)
        indices = [int(c.index) for c in markers.coordinate_identifier]
        if len(set(indices)) != len(indices):
            raise ValueError("Duplicate marker coordinate identifiers")
        missing = set(MARKER_COORDINATES) - set(indices)
        if missing:
            raise ValueError(
                f"Missing marker coordinate identifiers: {sorted(missing)}"
            )
        positions = array(markers.positions, "IDS marker positions", ndim=2)
        if positions.shape[1] != len(indices):
            raise ValueError("Marker coordinate count does not match positions")
        columns = [indices.index(i) for i in MARKER_COORDINATES]
        result.append(
            Markers(
                positions[:, columns],
                markers.weights,
                *_species(population.species),
                time=scalar(markers.time, "marker time"),
                rate=name == "distribution_sources",
            )
        )
    if not result:
        raise ValueError("IDS contains no marker populations")
    if len({p.time for p in result}) != 1:
        raise ValueError("All marker populations must use the same initial time")
    return result


def fill_species(species, markers):
    species.type.index = 2
    species.type.name = "ion"
    species.ion.element.resize(1)
    species.ion.element[0].a = markers.mass_number
    species.ion.element[0].z_n = markers.nuclear_charge
    species.ion.z_ion = markers.charge


def fill_marker_slice(target, positions, weights, time, *, diagnostics=None):
    target.time = time
    ncoord = 5 + (len(diagnostics) if diagnostics else 0)
    target.coordinate_identifier.resize(ncoord)
    for c, index, name in zip(
        target.coordinate_identifier, MARKER_COORDINATES, MARKER_NAMES
    ):
        c.index = index
        c.name = name
    values = np.asarray(positions).copy()
    if diagnostics:
        for k, (name, (index, description, data)) in enumerate(diagnostics.items(), 5):
            coord = target.coordinate_identifier[k]
            coord.index = index
            coord.name = name
            coord.description = description
            values = np.column_stack((values, data))
    target.positions = values
    target.weights = weights


def write_markers(populations, *, dd_version="4.1.1"):
    """Create a validated IDS without reinterpreting source/population weights."""
    populations = list(populations)
    if not populations or len({p.rate for p in populations}) != 1:
        raise ValueError("Provide populations with a common weight convention")
    times = {p.time for p in populations}
    if len(times) != 1:
        raise ValueError("Marker populations must have a common time")
    name = "distribution_sources" if populations[0].rate else "distributions"
    ids = factory(dd_version).new(name)
    initialize(ids, [populations[0].time])
    target = ids.source if populations[0].rate else ids.distribution
    target.resize(len(populations))
    for item, population in zip(target, populations):
        item.gyro_type = 2
        fill_species(item.species, population)
        item.markers.resize(1)
        fill_marker_slice(
            item.markers[0], population.positions, population.weights, population.time
        )
    ids.validate()
    return ids
