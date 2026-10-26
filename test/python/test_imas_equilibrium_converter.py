"""Exercise the IMAS equilibrium adapter without optional IMAS schemas."""
import importlib.metadata as metadata
import json
from pathlib import Path
import numpy as np
import pytest

libneo = pytest.importorskip("libneo")
netCDF4 = pytest.importorskip("netCDF4")
pytest.importorskip("scipy")
from imas_equilibrium_fixture import (model, normalized_radius,
                                      solovev_equilibrium, toroidal_flux)
from pysimple.imas.equilibrium import _geqdsk_input, equilibrium_to_chartmap


@pytest.fixture(scope="module", autouse=True)
def pinned_native_package():
    """Require the installed optional dependency, including its native ABI."""
    import _efit_to_boozer
    dist = metadata.distribution("libneo")
    assert dist.version == "0.0.4.dev0"
    direct = dist.read_text("direct_url.json")
    assert direct, "Install the pinned libneo VCS dependency, not an unrelated wheel"
    origin = json.loads(direct)
    assert origin["url"] == "https://github.com/itpplasma/libneo.git"
    assert origin["vcs_info"]["commit_id"] == "b60770c90b88bb3914c707ff82f0c161990ceffd"
    files = {Path(dist.locate_file(x)).resolve() for x in dist.files}
    assert Path(libneo.__file__).resolve() in files
    assert Path(_efit_to_boozer.__file__).resolve() in files


@pytest.mark.parametrize("version", ["3.42.2", "4.1.1"])
@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("field", [-1, 1])
def test_cocos_and_chartmap_extent(tmp_path, version, sign, field):
    eq = solovev_equilibrium(version, sign=sign, field=field)
    data = _geqdsk_input(eq)
    grid = eq.time_slice[0].profiles_2d[0].grid
    np.testing.assert_allclose(data["PsiVs"],
        sign*model()._psi(grid.dim1[None, :], grid.dim2[:, None]), atol=2e-15)
    np.testing.assert_allclose(np.sign(data["fprof"]), field, rtol=0, atol=0)
    path = equilibrium_to_chartmap(eq, tmp_path/"solovev.nc", nstep=360,
        nlabel=32, ntheta_int=64, nsurfmax=400, nsurf=64, nrho=24, ntheta=48)
    with netCDF4.Dataset(path) as ds:
        x, y, z = [np.asarray(ds.variables[n][:]) for n in "xyz"]
        flux = float(ds.torflux)
        assert ds.imas_cocos == (11 if version.startswith("3.") else 17)
    assert .98 < normalized_radius(np.hypot(x, y)*.01, z*.01).max() < 1.01
    np.testing.assert_allclose(flux, -toroidal_flux(field=field)*1e8, rtol=.01)


def test_rejects_original_non_equilibrium_toy(tmp_path):
    eq = solovev_equilibrium()
    t = eq.time_slice[0]
    r, z, s = np.linspace(8.7, 11.3, 65), np.linspace(-1.3, 1.3, 65), np.linspace(0, 1, 65)
    p2, p1, q = t.profiles_2d[0], t.profiles_1d, t.global_quantities
    p2.grid.dim1, p2.grid.dim2 = r, z
    p2.psi = np.pi*((r[:, None]-10)**2+z[None, :]**2)
    p1.psi, p1.f, p1.pressure = np.pi*s, np.full(65, 20.), np.zeros(65)
    p1.q = -20/np.sqrt(100-s)
    q.psi_axis, q.psi_boundary, q.magnetic_axis.r, q.magnetic_axis.z = 0., np.pi, 10., 0.
    q.ip = -2*np.pi/np.sqrt(99)/(4*np.pi*1e-7)
    theta = np.linspace(0, 2*np.pi, 129)
    t.boundary.outline.r, t.boundary.outline.z = 10+np.cos(theta), np.sin(theta)
    eq.vacuum_toroidal_field.r0, eq.vacuum_toroidal_field.b0 = 10., [2.]
    with pytest.raises(RuntimeError, match=r"(?i)Grad.Shafranov|force.balance"):
        equilibrium_to_chartmap(eq, tmp_path/"invalid.nc", nstep=360,
            nlabel=32, ntheta_int=64, nsurfmax=400, nsurf=64, nrho=24, ntheta=48)
