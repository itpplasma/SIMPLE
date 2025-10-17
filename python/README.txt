Clean SIMPLE Python API
=======================

The ``simple`` package provides a modern, batch-oriented interface to the
high-performance SIMPLE Fortran backend.  It mirrors the structure-of-arrays
memory layout used in Fortran, enabling zero-copy interaction while keeping the
public API concise.

Key modules
-----------

``simple.__init__``
    Public entry point exporting :class:`ParticleBatch`, :class:`BatchResults`,
    :class:`SurfaceSampler`, :class:`VolumeSampler`, and the :func:`trace_orbits`
    convenience function.

``simple.particles``
    Structure-of-arrays container that wraps particle phase-space coordinates.
    Provides helpers for sampling and for constructing batches from raw numpy
    arrays returned by the Fortran samplers.

``simple.results``
    Immutable view over the Fortran output arrays with convenience methods for
    confinement analysis.

``simple.samplers``
    Lightweight wrappers over the Fortran ``samplers`` module, exposing a
    Pythonic interface for surface/volume sampling as well as file-based
    particle loading.

Quick start
-----------

When working directly from the repository, expose the package via
``PYTHONPATH`` (or install it in editable mode) so that ``import simple``
resolves correctly:

.. code-block:: bash

    export PYTHONPATH=$PWD/python:$PYTHONPATH

If you prefer an installed package, use ``python -m pip install .`` which will
build the bindings and expose the same ``simple`` module.

The high-level :class:`simple.SimpleSession` takes care of VMEC loading and
keeps the Fortran globals in sync:

.. code-block:: python

    import simple

    session = simple.SimpleSession(simple.ensure_example_vmec())
    batch = session.sample_surface(1024, surface=0.4)
    results = session.trace(batch, tmax=0.2, integrator="symplectic_midpoint")

    confined_fraction = results.confined_mask().mean()
    print(f"Confined fraction: {confined_fraction:.2%}")

Testing
-------

See ``test/python`` for pytest-based validation of the API.  All tests load the
same VMEC equilibrium through the shared fixture in ``test/conftest.py``.
