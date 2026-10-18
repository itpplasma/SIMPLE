"""Heat loads on STL facets in file order; see DOC/wall-heat-loads.md."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import struct

import numpy as np


def _read_stl(path: Path) -> np.ndarray:
    """Read ASCII/binary STL without welding or reordering triangles."""
    raw = path.read_bytes()

    def ascii_triangles():
        # CGAL reads whitespace-separated tokens, independent of line breaks.
        tokens = iter(raw.decode('ascii').split())
        triangles = []
        for token in tokens:
            if token != 'facet':
                continue
            vertices = []
            for token in tokens:
                if token == 'endfacet':
                    break
                if token == 'vertex':
                    try:
                        vertices.append([float(next(tokens)) for _ in range(3)])
                    except StopIteration as exc:
                        raise ValueError('truncated STL vertex') from exc
                    if len(vertices) > 3:
                        raise ValueError('STL facets must have three vertices')
            else:
                raise ValueError('unterminated STL facet')
            if len(vertices) != 3:
                raise ValueError('STL facets must have three vertices')
            triangles.append(vertices)
        return np.asarray(triangles, dtype=float).reshape(-1, 3, 3)

    # Match CGAL: try ASCII for a 'solid' prefix, then binary on failure.
    prefix = raw.lstrip(b' \t\r\n\v\f')
    triangles = np.empty((0, 3, 3))
    if prefix[:5] == b'solid' and prefix[5:6] in (b' ', b'\r', b'\n'):
        try:
            triangles = ascii_triangles()
        except (UnicodeDecodeError, ValueError):
            pass
    if triangles.size == 0:
        count = struct.unpack_from('<I', raw, 80)[0] if len(raw) >= 84 else 0
        if count > 0 and len(raw) >= 84 + 50 * count:
            dtype = np.dtype([('normal', '<f4', 3), ('vertices', '<f4', (3, 3)),
                              ('attribute', '<u2')])
            triangles = np.frombuffer(raw, dtype=dtype, count=count,
                                      offset=84)['vertices'].copy()
        else:
            try:
                triangles = ascii_triangles()
            except (UnicodeDecodeError, ValueError) as exc:
                raise ValueError('invalid or truncated STL') from exc
    if triangles.size == 0 or not np.isfinite(triangles).all():
        raise ValueError('STL must contain finite triangles')
    return triangles.astype(float)


def _particle_array(ds, name: str) -> np.ndarray:
    var = ds.variables[name]
    if var.dimensions != ('particle',):
        raise ValueError(f'{name} must have dimension (particle,)')
    data = var[:]
    if np.ma.is_masked(data):
        raise ValueError(f'{name} contains missing data')
    return np.asarray(data)


@dataclass
class TriangleHeatLoad:
    """Unsmoothened facet powers and fluxes. Geometry is in m, power in W.

    Construct with ``from_netcdf``. Marker weights describe the *whole* source
    population and are normalized by their sum, including markers with no hit.
    ``zend`` momentum squared supplies the final energy fraction E/E0.
    """

    triangles_m: np.ndarray
    area_m2: np.ndarray
    hit_count: np.ndarray
    power_W: np.ndarray
    power_density_W_m2: np.ndarray
    power_density_per_MW: np.ndarray
    incidence_angle_rad: np.ndarray
    hit_power_W: np.ndarray
    total_alpha_power_W: float
    n_markers: int
    stl_file: str

    @classmethod
    def from_netcdf(cls, results_file: str | Path, stl_file: str | Path,
                    total_alpha_power_MW: float, *,
                    marker_weights: np.ndarray | None = None) -> TriangleHeatLoad:
        """Read hits, final energies and wall units from SIMPLE results.

        Use ``marker_weight(particle)`` if present, otherwise equal weights.
        Explicit ``marker_weights`` override that variable. STL identification
        is checked by basename and, when recorded, facet count. The caller
        must supply the original mesh with unchanged facet order.
        """
        import netCDF4 as nc

        if not np.isfinite(total_alpha_power_MW) or total_alpha_power_MW <= 0:
            raise ValueError('total_alpha_power_MW must be finite and positive')
        stl_file = Path(stl_file)
        with nc.Dataset(results_file) as ds:
            index_var = ds.variables['wall_hit_triangle']
            identity = index_var.getncattr('stl_file')
            if Path(identity).name != stl_file.name:
                raise ValueError('STL identification does not match results')
            if index_var.getncattr('start_index') != 0:
                raise ValueError('wall_hit_triangle must be 0-based')
            units = ds.getncattr('wall_units').strip().lower()
            scales = {'m': 1.0, 'mm': 0.001}
            if units not in scales:
                raise ValueError(f'unsupported wall_units: {units}')
            triangles = _read_stl(stl_file) * scales[units]
            if 'stl_facet_count' in index_var.ncattrs():
                if index_var.stl_facet_count != len(triangles):
                    raise ValueError('STL facet count does not match results')
            index = _particle_array(ds, 'wall_hit_triangle')
            flags = _particle_array(ds, 'wall_hit')
            if not np.isin(flags, [0, 1]).all():
                raise ValueError('wall_hit must be 0 or 1')
            hit = flags == 1
            n = len(index)
            if flags.shape != (n,) or n == 0:
                raise ValueError('wall hits must describe a nonempty population')
            if not np.issubdtype(index.dtype, np.integer):
                raise ValueError('wall_hit_triangle must contain integers')
            if (np.any(index[~hit] != -1) or np.any(index[hit] < 0)
                    or np.any(index[hit] >= len(triangles))):
                raise ValueError('wall_hit_triangle disagrees with wall_hit or mesh')
            zend = ds.variables['zend']
            if set(zend.dimensions) != {'phase', 'particle'}:
                raise ValueError('zend must have phase and particle dimensions')
            momentum = np.take(zend[:], 3, axis=zend.dimensions.index('phase'))
            if momentum.shape != (n,) or np.ma.is_masked(momentum[hit]):
                raise ValueError('invalid final momentum')
            # Empty masked-array reductions are masked, even without missing data.
            momentum = np.asarray(momentum)
            if not np.isfinite(momentum[hit]).all():
                raise ValueError('invalid final momentum')
            if np.any(momentum[hit] < 0):
                raise ValueError('final momentum magnitude must be nonnegative')
            angles = _particle_array(ds, 'wall_hit_angle_rad')[hit]
            if (not np.isfinite(angles).all() or np.any(angles < 0)
                    or np.any(angles > np.pi / 2)):
                raise ValueError('incidence angles must be in [0, pi/2] radians')
            if marker_weights is None:
                marker_weights = (_particle_array(ds, 'marker_weight')
                                  if 'marker_weight' in ds.variables
                                  else np.ones(n))
            weights = np.asarray(marker_weights, dtype=float)
            if (weights.shape != (n,) or not np.isfinite(weights).all()
                    or np.any(weights < 0) or not np.isfinite(weights.sum())
                    or weights.sum() <= 0):
                raise ValueError('marker weights must be finite, nonnegative, '
                                 'one per marker, with positive sum')

        area = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                       triangles[:, 2] - triangles[:, 0]), axis=1) / 2
        if not np.isfinite(area).all():
            raise ValueError('STL contains nonfinite triangle areas')
        if np.any(area[index[hit]] == 0):
            raise ValueError('recorded hit on a degenerate STL triangle')
        fractions = weights[hit] / weights.sum() * momentum[hit] ** 2
        total_W = float(total_alpha_power_MW) * 1e6
        count = np.bincount(index[hit], minlength=len(area))
        fraction = np.bincount(index[hit], weights=fractions, minlength=len(area))
        flux_per_MW = np.divide(fraction * 1e6, area,
                               out=np.zeros_like(area), where=area > 0)
        return cls(triangles, area, count, fraction * total_W,
                   flux_per_MW * float(total_alpha_power_MW), flux_per_MW,
                   angles, fractions * total_W, total_W, n, stl_file.name)

    @property
    def n_hits(self) -> int:
        return int(self.hit_count.sum())

    @property
    def deposited_power_W(self) -> float:
        return float(self.power_W.sum())

    @property
    def peak_W_m2(self) -> float:
        return float(self.power_density_W_m2.max())

    @property
    def p95_W_m2(self) -> float:
        """Area-weighted 95th percentile over all facets, including unhit area."""
        order = np.argsort(self.power_density_W_m2)
        cumulative = np.cumsum(self.area_m2[order])
        i = np.searchsorted(cumulative, 0.95 * cumulative[-1])
        return float(self.power_density_W_m2[order[i]])

    def wetted_area_m2(self, threshold_W_m2: float = 0.0) -> float:
        """Wall area with flux strictly above a nonnegative threshold."""
        if not np.isfinite(threshold_W_m2) or threshold_W_m2 < 0:
            raise ValueError('threshold must be finite and nonnegative')
        return float(self.area_m2[self.power_density_W_m2 > threshold_W_m2].sum())

    def incidence_distribution(self, bins: int = 18) -> dict[str, np.ndarray]:
        """Count and deposited-power histograms; angle 0 is normal incidence."""
        count, edges = np.histogram(self.incidence_angle_rad, bins=bins,
                                    range=(0, np.pi / 2))
        power, _ = np.histogram(self.incidence_angle_rad, bins=edges,
                                weights=self.hit_power_W)
        return {'edges_rad': edges, 'hit_count': count, 'power_W': power}

    def summary(self, threshold_W_m2: float = 0.0) -> dict:
        """JSON-serializable metrics, with units in their names."""
        distribution = self.incidence_distribution()
        return {'stl_file': self.stl_file, 'n_triangles': len(self.area_m2),
                'n_markers': self.n_markers, 'n_hits': self.n_hits,
                'total_alpha_power_W': self.total_alpha_power_W,
                'deposited_power_W': self.deposited_power_W,
                'peak_W_m2': self.peak_W_m2, 'p95_W_m2': self.p95_W_m2,
                'threshold_W_m2': threshold_W_m2,
                'wetted_area_m2': self.wetted_area_m2(threshold_W_m2),
                'incidence_distribution': {k: v.tolist()
                                           for k, v in distribution.items()}}

    def _cell_data(self) -> dict[str, np.ndarray]:
        return {'triangle_index': np.arange(len(self.area_m2)),
                'area_m2': self.area_m2, 'hit_count': self.hit_count,
                'power_W': self.power_W,
                'power_density_W_m2': self.power_density_W_m2,
                'power_density_W_m2_per_MW': self.power_density_per_MW}

    def to_vtk(self, path: str | Path) -> Path:
        """Write legacy ASCII VTK POLYDATA with facet cell data for ParaView."""
        path = Path(path)
        n = len(self.area_m2)
        with path.open('w') as out:
            out.write('# vtk DataFile Version 3.0\nSIMPLE wall heat load (m, W)\n'
                      f'ASCII\nDATASET POLYDATA\nPOINTS {3*n} double\n')
            np.savetxt(out, self.triangles_m.reshape(-1, 3), fmt='%.17g')
            out.write(f'POLYGONS {n} {4*n}\n')
            np.savetxt(out, np.column_stack((np.full(n, 3),
                                            np.arange(3*n).reshape(n, 3))), fmt='%d')
            out.write(f'CELL_DATA {n}\n')
            for name, values in self._cell_data().items():
                out.write(f'SCALARS {name} double 1\nLOOKUP_TABLE default\n')
                np.savetxt(out, values, fmt='%.17g')
        return path

    def to_netcdf(self, path: str | Path) -> Path:
        """Write triangle geometry and cell data, indexed as in the source STL."""
        import netCDF4 as nc

        path = Path(path)
        with nc.Dataset(path, 'w') as ds:
            ds.stl_file = self.stl_file
            ds.total_alpha_power_W = self.total_alpha_power_W
            for name, size in [('triangle', len(self.area_m2)), ('vertex', 3),
                               ('xyz', 3)]:
                ds.createDimension(name, size)
            vertices = ds.createVariable('vertices', 'f8', ('triangle', 'vertex', 'xyz'))
            vertices.units = 'm'
            vertices[:] = self.triangles_m
            units = ['1', 'm2', '1', 'W', 'W m-2', 'W m-2 MW-1']
            for (name, values), unit in zip(self._cell_data().items(), units):
                var = ds.createVariable(name, 'i8' if values.dtype.kind in 'iu'
                                        else 'f8', ('triangle',))
                var.units = unit
                var[:] = values
        return path

    def to_stl(self, path: str | Path) -> tuple[Path, Path]:
        """Write binary STL (m) plus ``.stl.nc`` cell-data sidecar.

        STL has no standard cell-data storage. The sidecar preserves exact
        indexed geometry and scalar fields for downstream engineering tools.
        """
        path = Path(path)
        normals = np.cross(self.triangles_m[:, 1] - self.triangles_m[:, 0],
                           self.triangles_m[:, 2] - self.triangles_m[:, 0])
        np.divide(normals, 2 * self.area_m2[:, None], out=normals,
                  where=self.area_m2[:, None] > 0)
        with path.open('wb') as out:
            out.write(b'SIMPLE wall heat load; coordinates in m'.ljust(80, b'\0'))
            out.write(struct.pack('<I', len(self.area_m2)))
            for normal, vertices in zip(normals, self.triangles_m):
                out.write(struct.pack('<12fH', *normal, *vertices.ravel(), 0))
        sidecar = self.to_netcdf(path.with_suffix(path.suffix + '.nc'))
        return path, sidecar
