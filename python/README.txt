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

.. code-block:: python

    import simple

    vmec = "wout.nc"
    simple.load_vmec(vmec)

    sampler = simple.SurfaceSampler(vmec)
    batch = simple.ParticleBatch.from_fortran_arrays(
        sampler.sample_surface_fieldline(1024, s=0.4)
    )

    results = simple.trace_orbits(
        batch,
        tmax=0.2,
        integrator="symplectic_midpoint",
    )

    confined_fraction = results.confined_mask().mean()
    print(f"Confined fraction: {confined_fraction:.2%}")

Testing
-------

See ``test/python`` for pytest-based validation of the API.  All tests load the
same VMEC equilibrium through the shared fixture in ``test/conftest.py``.
