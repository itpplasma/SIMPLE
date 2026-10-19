"""Bulk-plasma IDS profiles for SIMPLE's two-ion collision operator."""

import numpy as np

from ._common import array, check_ids, scalar, select_slice


def plasma_profiles(ids, *, time=None):
    """Return ``(table, config)`` in SI/eV without fitting the IDS profiles.

    SIMPLE derives electron density from quasineutrality. Reject inconsistent
    density, rotation, mixed charge states and more than two bulk species.
    The native collision operator samples normalized toroidal flux ``s`` for
    both VMEC and Boozer chartmaps; IDS rho_tor_norm is squared explicitly.
    """
    check_ids(ids, "core_profiles")
    p = select_slice(ids.profiles_1d, time)
    rho = array(p.grid.rho_tor_norm, "core_profiles rho_tor_norm")
    if len(rho) < 2 or rho[0] != 0 or rho[-1] != 1 or np.any(np.diff(rho) <= 0):
        raise ValueError("core_profiles rho_tor_norm must increase from 0 to 1")
    if not 1 <= len(p.ion) <= 2:
        raise ValueError("SIMPLE collisions require one or two bulk ion species")

    def profile(node, name):
        values = array(node, name, nonnegative=True)
        if len(values) != len(rho):
            raise ValueError(f"{name} length differs from rho_tor_norm")
        return values

    table = np.zeros((6, len(rho)), order="F")
    table[0] = rho**2
    table[1] = profile(p.electrons.temperature, "electron temperature")
    config = {"swcoll": True, "am2": 1.0, "Z2": 1.0}
    ne = np.zeros_like(rho)
    for i, ion in enumerate(p.ion):
        if len(ion.element) != 1 or len(ion.state):
            raise ValueError("Bulk ions must have one element and one charge state")
        mass = scalar(ion.element[0].a, "bulk ion mass number")
        charge = scalar(ion.z_ion, "bulk ion charge")
        if mass <= 0 or charge <= 0:
            raise ValueError("Bulk ion mass and charge must be positive")
        for rotation in (ion.rotation_frequency_tor,):
            if rotation.has_value and np.any(np.asarray(rotation) != 0):
                raise ValueError("Rotating plasma is unsupported by SIMPLE collisions")
        if ion.density_fast.has_value and np.any(np.asarray(ion.density_fast) != 0):
            raise ValueError("Bulk plasma must not include a fast-ion density")
        density = profile(
            ion.density_thermal if ion.density_thermal.has_value else ion.density,
            "bulk ion density",
        )
        table[2 + i] = profile(ion.temperature, "bulk ion temperature")
        table[4 + i] = density
        ne += charge * density
        config.update({f"am{i + 1}": mass, f"Z{i + 1}": charge})
    if len(p.ion) == 1:
        # Native collision coefficients still evaluate both temperatures.
        table[3] = table[2]
    electrons = profile(p.electrons.density, "electron density")
    if not np.allclose(ne, electrons, rtol=1e-6, atol=0):
        raise ValueError("SIMPLE requires quasineutral electron/ion densities")
    if np.any(table[1:4] <= 0) or np.any(ne <= 0):
        raise ValueError("Collision profiles require positive temperatures/density")
    return table, config
