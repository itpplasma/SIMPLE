Overview
========

The SIMPLE Python package provides a thin, batch-oriented wrapper around the
Fortran core.  The public API is intentionally small to reflect the underlying
high-performance code path.

Key concepts
------------

* :class:`simple.ParticleBatch` is a structure-of-arrays container that holds
  particle phase-space coordinates and interoperates directly with NumPy.
* :class:`simple.SimpleSession` encapsulates VMEC loading together with the
  most common operations (sampling, tracing, classification) so you don't need
  to manage global state manually.
* :class:`simple.BatchResults` exposes immutable access to simulation outputs,
  including convenience helpers for confinement analysis.
* :func:`simple.trace_orbits` runs the Fortran integrator for a prepared batch.
* :func:`simple.classify_fast` and
  :meth:`simple.SimpleSession.classify_fast` run the fast classifiers and, if
  requested, emit the legacy ``fort.*`` diagnostic files in a controlled
  directory.
* :class:`simple.SurfaceSampler` and :class:`simple.VolumeSampler` provide
  access to the Fortran sampling routines while keeping data in the new API
  types.

The documentation tree links to the API reference and to live examples that are
rendered directly from the scripts in :mod:`examples/`.  The recommended
workflow is:

1. Ensure a VMEC file is available via :func:`simple.ensure_example_vmec` or by
   supplying your own.
2. Create a :class:`simple.SimpleSession` with that VMEC file.
3. Use the session to sample particles, trace orbits, and perform
   classification, optionally writing the fort-style diagnostics for
   cross-checks with ``simple.x``.
