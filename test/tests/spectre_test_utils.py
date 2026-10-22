"""Extract the unchanged axis volume for single-volume orbit tests."""

from contextlib import contextmanager
from pathlib import Path
import tempfile

import numpy as np
from netCDF4 import Dataset


@contextmanager
def axis_volume(source):
    """Keep volume 1's field coefficients and its two bounding surfaces exactly."""
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / "axis_volume.h5"
        with Dataset(source, "r") as src, Dataset(target, "w") as dst:
            src.set_auto_mask(False)
            groups = {
                "input/physics": dst.createGroup("input").createGroup("physics"),
                "output": dst.createGroup("output"),
                "vector_potential": dst.createGroup("vector_potential"),
            }

            def write(key, values):
                values = np.asarray(values)
                parent, name = key.rsplit("/", 1)
                group = groups[parent]
                dims = tuple(f"{name}_d{i}" for i in range(values.ndim))
                for dim, size in zip(dims, values.shape):
                    group.createDimension(dim, size)
                group.createVariable(name, values.dtype, dims)[:] = values

            for key in ("input/physics/Igeometry", "input/physics/Nfp",
                        "input/physics/Mpol",
                        "input/physics/Ntor", "input/physics/phiedge",
                        "output/mn", "output/im", "output/in"):
                # phiedge is SPECTRE's global edge toroidal flux, kept unchanged
                # because real SPECTRE output always carries it and libneo's
                # reader requires it.
                write(key, src[key][:])
            for key in ("input/physics/Nvol", "output/Mvol"):
                write(key, [1])
            for key in ("input/physics/Lrad", "input/physics/pressure",
                        "input/physics/tflux", "input/physics/pflux", "output/mu"):
                write(key, src[key][:1])
            for key in ("Rbc", "Rbs", "Zbc", "Zbs"):
                write("output/" + key, src["output/" + key][:2])
            nrad = int(src["input/physics/Lrad"][0]) + 1
            for key in ("Ate", "Ato", "Aze", "Azo"):
                write("vector_potential/" + key,
                      src["vector_potential/" + key][:, :nrad])
        yield str(target)
