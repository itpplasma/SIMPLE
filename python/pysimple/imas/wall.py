"""IMAS Cartesian triangle walls in the same convention as ASCOT5."""

from pathlib import Path
import numpy as np

from ._common import array, check_ids, select_slice


def wall_to_stl(
    ids, output, *, time=None, description=0, space_index=0, toroidal_points=128
):
    """Write an ASCII STL in metres from one wall.description_ggd grid.

    IMAS node indices are one based. Cartesian coordinate identifiers may be
    reordered; triangle geometry is preserved. SIMPLE must be built with CGAL.
    """
    check_ids(ids, "wall")
    if not len(ids.description_ggd):
        units = ids.description_2d[description].limiter.unit
        if len(units) != 1:
            raise ValueError("Provide one closed axisymmetric limiter outline")
        outline = units[0].outline
        r, z = array(outline.r, "limiter R"), array(outline.z, "limiter Z")
        if len(r) != len(z) or len(r) < 4 or np.any(r <= 0):
            raise ValueError("Limiter outline needs >=4 positive-R points")
        if not np.allclose([r[0], z[0]], [r[-1], z[-1]], rtol=0, atol=1e-10):
            raise ValueError("Limiter outline must be explicitly closed")
        if not isinstance(toroidal_points, int) or toroidal_points < 8:
            raise ValueError("toroidal_points must be an integer >=8")
        r, z = r[:-1], z[:-1]
        if np.sum(r * np.roll(z, -1) - z * np.roll(r, -1)) < 0:
            r, z = r[::-1], z[::-1]
        phi = np.linspace(0, 2 * np.pi, toroidal_points, endpoint=False)
        nodes = np.stack(
            np.broadcast_arrays(
                r[None, :] * np.cos(phi[:, None]),
                r[None, :] * np.sin(phi[:, None]),
                z[None, :],
            ),
            axis=-1,
        )
        next_phi = np.roll(nodes, -1, axis=0)
        next_both = np.roll(next_phi, -1, axis=1)
        next_outline = np.roll(nodes, -1, axis=1)
        triangles = np.concatenate(
            (
                np.stack((nodes, next_phi, next_both), axis=-2).reshape(-1, 3, 3),
                np.stack((nodes, next_both, next_outline), axis=-2).reshape(-1, 3, 3),
            )
        )
        return _write_stl(triangles, output)
    static = int(ids.ids_properties.homogeneous_time) == 2
    grid = select_slice(
        ids.description_ggd[description].grid_ggd, None if static else time
    )
    space = grid.space[space_index]
    coordinates = [int(c.index) for c in space.coordinates_type]
    if sorted(coordinates) != [1, 2, 3]:
        raise ValueError("Wall coordinates must be Cartesian X/Y/Z [m] (IDs 1,2,3)")
    dimensions = space.objects_per_dimension
    if len(dimensions) < 3 or not len(dimensions[0].object):
        raise ValueError("Wall must contain nodes and triangle faces")
    nodes = np.stack([array(n.geometry, "wall node") for n in dimensions[0].object])
    if nodes.shape[1] != 3:
        raise ValueError("Wall nodes must have three coordinates")
    nodes = nodes[:, [coordinates.index(c) for c in [1, 2, 3]]]
    triangles = []
    for face in dimensions[2].object:
        indices = np.asarray(face.nodes, dtype=int)
        if indices.shape != (3,) or np.any(indices < 1) or np.any(indices > len(nodes)):
            raise ValueError(
                "Wall triangles require three valid one-based node indices"
            )
        triangles.append(nodes[indices - 1])
    if not triangles:
        raise ValueError("Wall contains no triangle faces")
    return _write_stl(triangles, output)


def _write_stl(triangles, output):
    triangles = np.asarray(triangles)
    normals = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    norms = np.linalg.norm(normals, axis=1)
    if np.any(norms == 0):
        raise ValueError("Wall contains a degenerate triangle")
    with Path(output).open("w") as stream:
        stream.write("solid imas_wall\n")
        for triangle, normal, norm in zip(triangles, normals, norms):
            normal = normal / norm
            stream.write(
                "facet normal " + " ".join(map(str, normal)) + "\nouter loop\n"
            )
            for node in triangle:
                stream.write("vertex " + " ".join(map(str, node)) + "\n")
            stream.write("endloop\nendfacet\n")
        stream.write("endsolid imas_wall\n")
    return Path(output)
