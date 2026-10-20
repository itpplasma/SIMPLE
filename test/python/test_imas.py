"""Real IMAS schemas, physical units and native actor integration."""

import numpy as np
import pytest

imas = pytest.importorskip("imas")

from pysimple.imas import (Markers, plasma_profiles, read_markers, run,
                          run_entry, write_markers)
from pysimple.imas._common import initialize
from pysimple.imas.equilibrium import _geqdsk_input, equilibrium_to_chartmap


def circular_equilibrium(version="4.1.1", *, sign=1, field=1, aspect=10):
    """Divergence-free circular flux surfaces, independent of native geometry.

    In physical cylindrical coordinates psi=C*((R-R0)**2+Z**2), F=R0*B0.
    This analytic field need not satisfy scalar-pressure force balance.
    """
    eq = imas.IDSFactory(version=version).new("equilibrium")
    initialize(eq, [0])
    eq.vacuum_toroidal_field.r0 = aspect
    eq.vacuum_toroidal_field.b0 = [2*field]
    eq.time_slice.resize(1)
    t = eq.time_slice[0]
    t.time = 0
    scale = (-1 if version.startswith("3.") else 1)*2*np.pi
    r = np.linspace(aspect-1.3, aspect+1.3, 65)
    z = np.linspace(-1.3, 1.3, 65)
    c = sign*.5
    t.profiles_2d.resize(1)
    p2 = t.profiles_2d[0]
    p2.grid_type.index = 1
    p2.grid.dim1, p2.grid.dim2 = r, z
    p2.psi = scale*c*((r[:, None]-aspect)**2+z[None, :]**2)
    setattr(p2, "b_field_tor" if version.startswith("3.") else "b_field_phi",
            field*2*aspect/r[:, None]*np.ones((1, len(z))))
    s = np.linspace(0, 1, len(r))
    p1 = t.profiles_1d
    p1.psi = scale*c*s
    p1.f = np.full(len(s), field*2*aspect)
    p1.pressure = np.zeros(len(s))
    p1.q = -field*2*aspect/(2*c*np.sqrt(aspect**2-s))
    p1.rho_tor_norm = np.sqrt((aspect-np.sqrt(aspect**2-s)) /
                             (aspect-np.sqrt(aspect**2-1)))
    q = t.global_quantities
    q.psi_axis = 0
    q.psi_boundary = scale*c
    q.magnetic_axis.r, q.magnetic_axis.z = aspect, 0
    q.ip = -4*np.pi*c/np.sqrt(aspect**2-1)/(4*np.pi*1e-7)
    theta = np.linspace(0, 2*np.pi, 129)
    t.boundary.outline.r = aspect+np.cos(theta)
    t.boundary.outline.z = np.sin(theta)
    eq.validate()
    return eq


def source(*, rate=False, version="4.1.1"):
    p = Markers(np.array([[10.5, .3, 0, 1e4, .8],
                          [10.6, -.7, 0, 2e4, -.6]]),
                np.array([2., 3.]), 2, 1, 1, rate=rate)
    return write_markers([p], dd_version=version)


@pytest.mark.parametrize("version", ["3.42.2", "4.1.1"])
@pytest.mark.parametrize("rate", [True, False])
def test_markers_real_schema_and_netcdf(tmp_path, version, rate):
    ids = source(rate=rate, version=version)
    path = tmp_path/"markers.nc"
    with imas.DBEntry(str(path), "w", dd_version=version) as db:
        db.put(ids)
    with imas.DBEntry(str(path), "r", dd_version=version) as db:
        loaded = db.get(ids.metadata.name)
    result = read_markers(loaded)[0]
    np.testing.assert_array_equal(result.positions, read_markers(ids)[0].positions)
    np.testing.assert_array_equal(result.weights, [2, 3])
    assert result.rate == rate


def test_markers_identifiers_and_validation():
    ids = source()
    m = ids.distribution[0].markers[0]
    order = [3, 1, 4, 0, 2]
    identifiers = [int(c.index) for c in m.coordinate_identifier]
    original = m.positions.value.copy()
    for c, i in zip(m.coordinate_identifier, order):
        c.index = identifiers[i]
    m.positions = original[:, order]
    np.testing.assert_array_equal(read_markers(ids)[0].positions, original)
    m.coordinate_identifier[0].index = 4
    with pytest.raises(ValueError, match="Duplicate"):
        read_markers(ids)
    with pytest.raises(ValueError, match="No IDS slice"):
        read_markers(source(), time=1)
    with pytest.raises(ValueError, match="source_interval"):
        run(source(rate=True), field="unused", duration=1e-6)


@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("version", ["3.42.2", "4.1.1"])
def test_cocos_conversion_against_analytic_field(sign, version):
    eq = circular_equilibrium(version, sign=sign)
    data = _geqdsk_input(eq)
    r, z = eq.time_slice[0].profiles_2d[0].grid.dim1, eq.time_slice[0].profiles_2d[0].grid.dim2
    expected = sign*.5*((np.asarray(r)[None, :]-10)**2+np.asarray(z)[:, None]**2)
    np.testing.assert_allclose(data["PsiVs"], expected, atol=2e-15)
    np.testing.assert_allclose(data["fprof"], 20)
    np.testing.assert_allclose(data["qprof"], sign*20/np.sqrt(100-np.linspace(0, 1, 65)))


def core_profiles():
    ids = imas.IDSFactory(version="4.1.1").new("core_profiles")
    initialize(ids, [0])
    ids.profiles_1d.resize(1)
    p = ids.profiles_1d[0]
    p.time = 0
    p.grid.rho_tor_norm = [0, .25, 1]
    p.electrons.density = [5e19]*3
    p.electrons.temperature = [1e4, 8e3, 1e3]
    p.ion.resize(1)
    ion = p.ion[0]
    ion.element.resize(1)
    ion.element[0].a, ion.element[0].z_n = 2, 1
    ion.z_ion = 1
    ion.density_thermal = [5e19]*3
    ion.temperature = [8e3, 6e3, 1e3]
    ids.validate()
    return ids


def test_bulk_profiles_no_fitting_or_species_loss():
    ids = core_profiles()
    table, config = plasma_profiles(ids)
    np.testing.assert_array_equal(table[0], [0, .0625, 1])
    np.testing.assert_array_equal(table[1], [1e4, 8e3, 1e3])
    assert config["am1"] == 2 and config["Z1"] == 1
    ids.profiles_1d[0].electrons.density = [1e19]*3
    with pytest.raises(ValueError, match="quasineutral"):
        plasma_profiles(ids)


@pytest.mark.parametrize("field", [1, -1])
@pytest.mark.parametrize("sign", [1, -1])
@pytest.mark.parametrize("version", ["3.42.2", "4.1.1"])
def test_chartmap_extent_all_flux_directions(tmp_path, version, sign, field):
    """libneo's scan stops at a signed flux; a reversed sign gives ~3 mm."""
    netCDF4 = pytest.importorskip("netCDF4")
    pytest.importorskip("_efit_to_boozer")
    path = equilibrium_to_chartmap(
        circular_equilibrium(version, sign=sign, field=field), tmp_path/"c.nc",
        nstep=360, nlabel=32, ntheta_int=64, nsurfmax=400, nsurf=64, nrho=24,
        ntheta=48)
    with netCDF4.Dataset(path) as ds:
        x, y, z = [np.asarray(ds.variables[n][:]) for n in "xyz"]
        torflux = float(ds.torflux)
    minor = np.sqrt((np.hypot(x, y)*.01-10)**2+(z*.01)**2)
    assert .98 < minor.max() < 1.01
    # Toroidal flux per radian of the unit disc at R*B_phi=20 T m, in G cm^2;
    # SIMPLE's chartmap convention gives it the sign opposite to B_phi.
    np.testing.assert_allclose(torflux, -field*20*(10-np.sqrt(99))*1e8, rtol=.01)


@pytest.fixture(scope="module", params=[1, -1])
def circular_chartmap(tmp_path_factory, request):
    pytest.importorskip("_efit_to_boozer")
    path = tmp_path_factory.mktemp("imas-field")/"circular.nc"
    equilibrium_to_chartmap(circular_equilibrium(sign=request.param), path, nstep=360,
        nlabel=32, ntheta_int=64, nsurfmax=400, nsurf=64, nrho=24, ntheta=48)
    return path


def test_native_physical_coordinate_roundtrip(circular_chartmap):
    import pysimple
    pysimple.init(circular_chartmap, ntestpart=2, trace_time=1e-6)
    cylindrical = np.array([[10.5, .3, 0], [10.3, -.7, .4]]).T
    reference = pysimple.cylindrical_to_reference(cylindrical)
    expected_rho = np.sqrt((10-np.sqrt(100-.25))/(10-np.sqrt(99)))
    np.testing.assert_allclose(reference[0], expected_rho, atol=5e-3, rtol=0)
    result = pysimple.reference_to_cylindrical(reference)
    np.testing.assert_allclose(result, cylindrical, atol=2e-7)
    with pytest.raises(ValueError, match="inverse failed"):
        pysimple.cylindrical_to_reference([-1, 0, 0])


def test_native_actor_particle_and_energy_totals(circular_chartmap, tmp_path):
    output = run(source(rate=True), field=circular_chartmap,
                 duration=1e-6, source_interval=.2, workdir=tmp_path)
    output.validate()
    d = output.distribution[0]
    np.testing.assert_allclose(d.markers[0].weights, [.4, .6])
    np.testing.assert_allclose(d.markers[1].weights, [.4, .6])
    np.testing.assert_allclose(d.markers[1].positions[:, 3], [1e4, 2e4], rtol=1e-6)
    assert float(d.global_quantities[1].particles_n) == pytest.approx(1)
    assert float(d.global_quantities[1].energy) == pytest.approx(1.6e4*1.602176634e-19, rel=1e-6)
    path = tmp_path/"result.nc"
    with imas.DBEntry(str(path), "w", dd_version="4.1.1") as db:
        db.put(output)
    with imas.DBEntry(str(path), "r", dd_version="4.1.1") as db:
        loaded = db.get("distributions")
    np.testing.assert_array_equal(loaded.distribution[0].markers[1].positions,
                                  d.markers[1].positions)


def test_actor_rejects_invalid_resolved_native_output(circular_chartmap, monkeypatch):
    import pysimple
    native_trace = pysimple.trace_parallel

    def corrupt_final_state(*args, **kwargs):
        result = native_trace(*args, **kwargs)
        # A native LCFS/completion code alone does not guarantee valid states.
        result["final_positions"][0, 0] = 2.0
        return result

    monkeypatch.setattr(pysimple, "trace_parallel", corrupt_final_state)
    with pytest.raises(RuntimeError, match="Invalid final states"):
        run(source(), field=circular_chartmap, duration=1e-6)


def test_database_actor_and_core_profiles(tmp_path):
    input_path, output_path = tmp_path/"input.nc", tmp_path/"output.nc"
    alpha = read_markers(source(rate=True))[0]
    alpha.mass_number, alpha.nuclear_charge, alpha.charge = 4, 2, 2
    with imas.DBEntry(str(input_path), "w", dd_version="4.1.1") as entry:
        for ids in [circular_equilibrium(), write_markers([alpha]), core_profiles()]:
            entry.put(ids)
    output = run_entry(input_path, output_path, core_profiles=True,
        duration=1e-7, source_interval=.2, deterministic=True, workdir=tmp_path,
        conversion=dict(nstep=360, nlabel=32, ntheta_int=64,
                        nsurfmax=400, nsurf=64, nrho=24, ntheta=48))
    output.validate()
    positions = output.distribution[0].markers[1].positions
    assert np.all(np.isfinite(positions))
    assert np.all(np.asarray(positions)[:, 3] > 0)
    assert output_path.is_file()
    # The actor's temporary field file is gone; the loaded physical map must
    # remain usable without reopening its metadata on every coordinate call.
    import pysimple
    assert np.all(np.isfinite(pysimple.reference_to_cylindrical([.5, 0, 0])))


def test_wall_ids_rejected_before_tracing(circular_chartmap, tmp_path):
    wall = imas.IDSFactory(version="4.1.1").new("wall")
    with pytest.raises(ValueError, match="not yet supported.*LCFS.*Omit the wall"):
        run(source(rate=True), field=circular_chartmap, wall=wall,
            duration=1e-6, source_interval=.2, workdir=tmp_path)
    with pytest.raises(ValueError, match="not yet supported"):
        run_entry(tmp_path/"in.nc", tmp_path/"out.nc", wall=True, duration=1e-6)
    assert not (tmp_path/"out.nc").exists()
