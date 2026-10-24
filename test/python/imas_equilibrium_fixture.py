"""Public analytic Solov'ev oracle from the pinned libneo test source.

Its polynomial obeys Delta*psi=(8B+2A)R²+2C=-mu0 R² p'-FF'.
The locator reads a fixture only and never changes dependency resolution.
"""
from functools import lru_cache
import importlib.util
import os
from pathlib import Path
from types import SimpleNamespace as NS
import numpy as np


@lru_cache(None)
def model():
    path = Path(os.environ["SIMPLE_TEST_LIBNEO_SOURCE_DIR"])/"python/tests/efit_solovev.py"
    spec = importlib.util.spec_from_file_location("simple_solovev_oracle", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def contour(level):
    m = model()
    angle = np.linspace(0, 2*np.pi, 1025)
    co, si = np.cos(angle), np.sin(angle)
    axis = m._axis()
    low, high = np.zeros_like(angle), np.full_like(angle, 1.5)
    for _ in range(65):
        rho = (low+high)/2
        outside = m._psi(axis[0]+rho*co, axis[1]+rho*si) > level
        high, low = np.where(outside, rho, high), np.where(outside, low, rho)
    return axis, (low+high)/2, co, si


def toroidal_flux(level=None, field=1):
    """Independent integral of F/R over the cross-section, per radian (SI)."""
    m = model()
    axis, edge, co, si = contour(m.PSI_EDGE if level is None else level)
    nodes, weights = np.polynomial.legendre.leggauss(48)
    rho = edge[:-1, None]*(nodes+1)/2
    r, z = axis[0]+rho*co[:-1, None], axis[1]+rho*si[:-1, None]
    f = field*np.sqrt((m.R0*m.B0)**2-4*m.C*(m._psi(r, z)-m.PSI_EDGE))
    return np.mean(np.sum(f/r*rho*weights*edge[:-1, None]/2, axis=1))


def normalized_radius(r, z):
    m = model()
    axis = m._psi(*m._axis())
    return np.sqrt(np.maximum(0, (m._psi(r, z)-axis)/(m.PSI_EDGE-axis)))


def expected_rho(r, z):
    levels = model()._psi(np.asarray(r), np.asarray(z))
    return np.sqrt(np.array([toroidal_flux(x) for x in levels.flat]).reshape(levels.shape)
                   / toroidal_flux())


def solovev_equilibrium(version="4.1.1", *, sign=1, field=1, schema=False):
    m = model()
    axis, radius, co, si = contour(m.PSI_EDGE)
    r, z = np.linspace(.9, 2.4, 65), np.linspace(-.9, .9, 65)
    flux = np.linspace(m._psi(*axis), m.PSI_EDGE, len(r))
    scale = (-1 if version.startswith("3.") else 1)*2*np.pi
    p2 = NS(grid_type=NS(index=1), grid=NS(dim1=r, dim2=z),
            psi=scale*sign*m._psi(r[:, None], z[None, :]))
    p1 = NS(psi=scale*sign*flux,
            f=field*np.sqrt((m.R0*m.B0)**2-4*m.C*(flux-m.PSI_EDGE)),
            pressure=-(8*m.B+2*m.A)/m.MU0*(flux-m.PSI_EDGE),
            q=np.array([-m._contour_oracle(x, sign, field)[0]
                        for x in flux[0]+np.maximum(np.linspace(0, 1, len(r)), 1e-6)
                        *(m.PSI_EDGE-flux[0])]))
    q = NS(psi_axis=scale*sign*flux[0], psi_boundary=scale*sign*m.PSI_EDGE,
           magnetic_axis=NS(r=axis[0], z=axis[1]),
           ip=m._contour_oracle(m.PSI_EDGE, sign, field)[1]*5e6)
    outline = NS(r=axis[0]+radius*co, z=axis[1]+radius*si)
    t = NS(time=0, profiles_2d=[p2], profiles_1d=p1, global_quantities=q,
           boundary=NS(outline=outline))
    vacuum = NS(r0=m.R0, b0=[field*m.B0])
    if not schema:
        return NS(_dd_version=version, metadata=NS(name="equilibrium"),
                  time_slice=[t], vacuum_toroidal_field=vacuum)
    import imas
    from pysimple.imas._common import initialize
    eq = imas.IDSFactory(version=version).new("equilibrium")
    initialize(eq, [0])
    eq.vacuum_toroidal_field.r0, eq.vacuum_toroidal_field.b0 = vacuum.r0, vacuum.b0
    eq.time_slice.resize(1)
    target = eq.time_slice[0]
    target.time = 0
    target.profiles_2d.resize(1)
    target.profiles_2d[0].grid_type.index = 1
    target.profiles_2d[0].grid.dim1, target.profiles_2d[0].grid.dim2 = r, z
    target.profiles_2d[0].psi = p2.psi
    for name in ("psi", "f", "pressure", "q"):
        setattr(target.profiles_1d, name, getattr(p1, name))
    for name in ("psi_axis", "psi_boundary", "ip"):
        setattr(target.global_quantities, name, getattr(q, name))
    target.global_quantities.magnetic_axis.r, target.global_quantities.magnetic_axis.z = axis
    target.boundary.outline.r, target.boundary.outline.z = outline.r, outline.z
    eq.validate()
    return eq
