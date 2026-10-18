"""Independent power/geometry oracles and the recorded CGAL tetrahedron fixture."""
import os
from pathlib import Path
import struct

import netCDF4 as nc
import numpy as np
import pytest

from pysimple.wall_heat_load import TriangleHeatLoad


# Unequal areas (0.5 and 2 m2) make flux and area quantiles distinguishable.
TRIANGLES = np.array([[[0, 0, 0], [1, 0, 0], [0, 1, 0]],
                      [[0, 0, 1], [2, 0, 1], [0, 2, 1]]], dtype=float)


def write_stl(path, triangles=TRIANGLES, binary=False):
    if binary:
        # A binary STL header may start with 'solid'; size determines the format.
        with path.open('wb') as f:
            f.write(b'solid binary fixture'.ljust(80, b'\0'))
            f.write(struct.pack('<I', len(triangles)))
            for tri in triangles:
                f.write(struct.pack('<12fH', 0, 0, 1, *tri.ravel(), 0))
    else:
        with path.open('w') as f:
            f.write('solid fixture\n')
            for tri in triangles:
                f.write('facet normal 0 0 1\nouter loop\n')
                for vertex in tri:
                    f.write('vertex ' + ' '.join(map(str, vertex)) + '\n')
                f.write('endloop\nendfacet\n')
            f.write('endsolid fixture\n')


def write_results(path, *, index=(1, 0, 1, -1), momentum=(1, 1, 1, 1),
                  weights=None, phase_first=False, units='m'):
    with nc.Dataset(path, 'w') as ds:
        ds.createDimension('particle', len(index))
        ds.createDimension('phase', 5)
        ds.wall_units = units
        v = ds.createVariable('wall_hit_triangle', 'i4', ('particle',))
        v.stl_file = 'wall.stl'
        v.start_index = 0
        v[:] = index
        ds.createVariable('wall_hit', 'i1', ('particle',))[:] = np.array(index) >= 0
        ds.createVariable('wall_hit_angle_rad', 'f8', ('particle',))[:] = [
            0, np.pi / 4, np.pi / 2, 0]
        shape = ('phase', 'particle') if phase_first else ('particle', 'phase')
        v = ds.createVariable('zend', 'f8', shape)
        z = np.zeros((5, len(index)))
        z[3] = momentum
        v[:] = z if phase_first else z.T
        if weights is not None:
            ds.createVariable('marker_weight', 'f8', ('particle',))[:] = weights


@pytest.fixture
def inputs(tmp_path):
    stl, results = tmp_path / 'wall.stl', tmp_path / 'results.nc'
    write_stl(stl)
    write_results(results)
    return results, stl


@pytest.mark.parametrize('binary', [False, True])
@pytest.mark.parametrize('units,scale', [('m', 1), ('cm', 100), ('mm', 1000)])
def test_conservation_and_known_bins(tmp_path, binary, units, scale):
    stl, results = tmp_path / 'wall.stl', tmp_path / 'results.nc'
    write_stl(stl, TRIANGLES * scale, binary=binary)
    write_results(results, index=[1, 0, 1, 0], units=units)
    load = TriangleHeatLoad.from_netcdf(results, stl, 2.0)
    np.testing.assert_allclose(load.area_m2, [0.5, 2], rtol=1e-14)
    np.testing.assert_array_equal(load.hit_count, [2, 2])
    np.testing.assert_allclose(load.power_W, [1e6, 1e6], rtol=1e-14)
    assert load.deposited_power_W == pytest.approx(2e6, rel=1e-14)
    np.testing.assert_allclose(load.power_density_W_m2, [2e6, 0.5e6])
    np.testing.assert_allclose(load.power_density_per_MW, [1e6, 0.25e6])


@pytest.mark.parametrize('phase_first', [False, True])
def test_weight_and_final_energy(inputs, phase_first):
    results, stl = inputs
    write_results(results, momentum=[0.5, 1, 0.25, 1],
                  weights=[1, 2, 3, 4], phase_first=phase_first)
    load = TriangleHeatLoad.from_netcdf(results, stl, 1)
    # All ten source weight units are in the denominator, even the confined four.
    np.testing.assert_allclose(load.power_W, [200000, 43750])
    assert load.n_hits == 3
    override = TriangleHeatLoad.from_netcdf(results, stl, 1,
                                           marker_weights=np.ones(4))
    np.testing.assert_allclose(override.power_W, [250000, 78125])


def test_metrics_and_incidence(inputs):
    load = TriangleHeatLoad.from_netcdf(*inputs, 1)
    assert load.peak_W_m2 == 500000
    assert load.p95_W_m2 == 500000
    assert load.wetted_area_m2() == 2.5
    assert load.wetted_area_m2(250000) == 0.5
    distribution = load.incidence_distribution(bins=3)
    np.testing.assert_array_equal(distribution['hit_count'], [1, 1, 1])
    np.testing.assert_allclose(distribution['power_W'], [250000] * 3)
    assert load.summary()['n_hits'] == 3
    # 80% of the area has lower flux, so the area p95 is the high-flux cell.
    assert load.p95_W_m2 != np.percentile(load.power_density_W_m2, 95)


def test_area_weighted_percentile(inputs):
    results, stl = inputs
    triangles = TRIANGLES.copy()
    triangles[0] *= 0.1  # high-flux area is less than 5% of the mesh
    write_stl(stl, triangles)
    load = TriangleHeatLoad.from_netcdf(results, stl, 1)
    assert load.peak_W_m2 == pytest.approx(50e6)
    assert load.p95_W_m2 == 250000


def test_no_hits(inputs):
    results, stl = inputs
    write_results(results, index=[-1] * 4)
    load = TriangleHeatLoad.from_netcdf(results, stl, 1)
    assert load.n_hits == 0
    assert load.deposited_power_W == load.peak_W_m2 == load.p95_W_m2 == 0
    assert load.wetted_area_m2() == 0
    assert load.incidence_distribution()['power_W'].sum() == 0


@pytest.mark.parametrize('weights', [[0, 0, 0, 0], [1, -1, 1, 1],
                                    [1, np.nan, 1, 1], [1, 2]])
def test_invalid_weights(inputs, weights):
    with pytest.raises(ValueError, match='marker weights'):
        TriangleHeatLoad.from_netcdf(*inputs, 1, marker_weights=weights)


@pytest.mark.parametrize('index,flag', [(8, 1), (-1, 1), (0, 0)])
def test_inconsistent_index(inputs, index, flag):
    results, stl = inputs
    with nc.Dataset(results, 'a') as ds:
        ds['wall_hit_triangle'][0] = index
        ds['wall_hit'][0] = flag
    with pytest.raises(ValueError, match='disagrees'):
        TriangleHeatLoad.from_netcdf(results, stl, 1)


def test_mesh_identity_and_degenerate_area(inputs):
    results, stl = inputs
    with nc.Dataset(results, 'a') as ds:
        ds['wall_hit_triangle'].stl_file = 'other.stl'
    with pytest.raises(ValueError, match='identification'):
        TriangleHeatLoad.from_netcdf(results, stl, 1)
    write_results(results)
    write_stl(stl, np.zeros_like(TRIANGLES))
    with pytest.raises(ValueError, match='degenerate'):
        TriangleHeatLoad.from_netcdf(results, stl, 1)


@pytest.mark.parametrize('power', [0, -1, np.inf, np.nan])
def test_invalid_power(inputs, power):
    with pytest.raises(ValueError, match='positive'):
        TriangleHeatLoad.from_netcdf(*inputs, power)


def test_exports(inputs, tmp_path):
    load = TriangleHeatLoad.from_netcdf(*inputs, 1)
    vtk_path = load.to_vtk(tmp_path / 'heat.vtk')
    # Read through VTK's public reader as an independent cell-data oracle.
    vtk = pytest.importorskip('vtk')
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(str(vtk_path))
    reader.ReadAllScalarsOn()
    reader.Update()
    mesh = reader.GetOutput()
    assert mesh.GetNumberOfCells() == 2
    for name, expected in {
        'triangle_index': [0, 1], 'area_m2': [0.5, 2], 'hit_count': [1, 2],
        'power_W': [250000, 500000],
        'power_density_W_m2': [500000, 250000],
        'power_density_W_m2_per_MW': [500000, 250000],
    }.items():
        actual = mesh.GetCellData().GetArray(name)
        np.testing.assert_allclose([actual.GetTuple1(i) for i in range(2)], expected)
    stl_path, sidecar = load.to_stl(tmp_path / 'heat.stl')
    stl_reader = vtk.vtkSTLReader()
    stl_reader.SetFileName(str(stl_path))
    stl_reader.Update()
    assert stl_reader.GetOutput().GetNumberOfCells() == 2
    with nc.Dataset(sidecar) as ds:
        np.testing.assert_allclose(ds['vertices'][:], TRIANGLES)
        np.testing.assert_allclose(ds['power_density_W_m2'][:], [500000, 250000])
        np.testing.assert_array_equal(ds['triangle_index'][:], [0, 1])


def test_recorded_wall_fixture():
    fixture = os.environ.get('SIMPLE_WALL_FIXTURE')
    if fixture is None:
        pytest.skip('requires CGAL wall fixture; run with make test and CGAL enabled')
    root = Path(fixture)
    load = TriangleHeatLoad.from_netcdf(root / 'tetra_wall_results.nc',
                                       root / 'tetra_m.stl', 1)
    with nc.Dataset(root / 'tetra_wall_results.nc') as ds:
        assert load.n_hits == np.count_nonzero(ds['wall_hit'][:]) == 4
        np.testing.assert_array_equal(ds['wall_hit_triangle'][:], [0, 1, 2, 3, -1])
        # Independent geometry oracle: each recorded point is in its STL plane.
        assert ds['wall_hit_cart'].units == 'cm'
        points = ds['wall_hit_cart'][:4] * 0.01
        normals = ds['wall_hit_normal_cart'][:4]
    residual = np.einsum('ij,ij->i', points - load.triangles_m[:, 0], normals)
    np.testing.assert_allclose(residual, 0, atol=1e-12)
    np.testing.assert_array_equal(load.hit_count, [1, 1, 1, 1])
    assert load.deposited_power_W == pytest.approx(800000, rel=1e-14)
