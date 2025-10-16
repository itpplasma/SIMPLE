Overview
========

The SIMPLE Python package provides a thin, batch-oriented wrapper around the
Fortran core.  The public API is intentionally small to reflect the underlying
high-performance code path.

Key concepts
------------

* :class:`simple.ParticleBatch` is a structure-of-arrays container that holds
  particle phase-space coordinates and interoperates directly with NumPy.
* :class:`simple.BatchResults` exposes immutable access to simulation outputs,
  including convenience helpers for confinement analysis.
* :func:`simple.trace_orbits` runs the Fortran integrator for a prepared batch.
* :func:`simple.classify_fast` runs the fast classifiers and returns their labels
  without generating legacy ``fort.*`` files.
* :class:`simple.SurfaceSampler` and :class:`simple.VolumeSampler` provide
  access to the Fortran sampling routines while keeping data in the new API
  types.

The documentation tree links to the API reference and to live examples that are
rendered directly from the scripts in :mod:`examples/`.
