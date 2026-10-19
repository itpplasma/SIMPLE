"""Shared IMAS schema checks and provenance for SIMPLE's actor."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
import numpy as np

DD_VERSION = "4.1.1"


def require_imas():
    try:
        import imas
    except ImportError as exc:
        raise ImportError("Install pysimple[imas] to use the IMAS interface") from exc
    return imas


def factory(dd_version=DD_VERSION):
    return require_imas().IDSFactory(version=dd_version)


def ids_version(ids):
    version = ids._dd_version
    if not (version.startswith("3.") or version.startswith("4.")):
        raise ValueError(f"Unsupported IMAS Data Dictionary version {version}")
    return version


def check_ids(ids, name):
    if ids.metadata.name != name:
        raise ValueError(f"Expected {name} IDS, got {ids.metadata.name}")
    ids_version(ids)


def array(node, name, *, ndim=1, nonnegative=False):
    values = np.asarray(node, dtype=float)
    if values.ndim != ndim or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must be a nonempty finite {ndim}-D array")
    if nonnegative and np.any(values < 0):
        raise ValueError(f"{name} must be nonnegative")
    return values.copy()


def scalar(node, name):
    value = float(node)
    if not np.isfinite(value) or value < -1e30:
        raise ValueError(f"{name} is missing or nonfinite")
    return value


def select_slice(slices, time=None):
    if not len(slices):
        raise ValueError("IDS contains no time slices")
    if time is None:
        return slices[0]
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    for item in slices:
        if np.isclose(float(item.time), time, rtol=0, atol=1e-10):
            return item
    raise ValueError(f"No IDS slice at time {time}; select/interpolate a slice first")


def initialize(ids, times):
    ids.ids_properties.homogeneous_time = 1
    ids.time = np.asarray(times, dtype=float)
    ids.code.name = "SIMPLE"
    try:
        ids.code.version = version("pysimple")
    except PackageNotFoundError:
        ids.code.version = "source checkout"
    ids.code.repository = "https://github.com/itpplasma/SIMPLE"
    ids.code.output_flag = np.ones(len(times), dtype=np.int32)
