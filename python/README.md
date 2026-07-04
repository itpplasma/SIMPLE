Clean SIMPLE Python API
=======================

The supported Python package is ``pysimple``. It wraps the Fortran backend built
by CMake/f90wrap and provides a module-level API for initialization, sampling,
tracing, and classification.

Environment setup
-----------------

For a fresh checkout, create the recommended local virtual environment:

.. code-block:: bash

    ./setup-venv.sh

Later, reactivate it with:

.. code-block:: bash

    source .venv/bin/activate

The editable install performed by ``setup-venv.sh`` exposes ``import pysimple``
directly from the repository checkout.

If you are working from the sibling benchmark checkout layout
(``../benchmark-simple-potato``, ``../SIMPLE``, ``../NEO-RT``), use the shared
environment in ``../benchmark-simple-potato/.venv`` instead of maintaining a
second venv here.

Quick start
-----------

First verify the bindings import cleanly:

.. code-block:: bash

    python -c "import pysimple; print('pysimple ok')"

.. code-block:: python

    import pysimple

    pysimple.init("wout.nc", deterministic=True, trace_time=5e-5, ntestpart=32)
    particles = pysimple.sample_surface(32, s=0.3)
    results = pysimple.trace_parallel(particles, integrator="midpoint")

    n_lost = (results["loss_times"] < 5e-5).sum()
    print(f"Lost particles: {n_lost}")

Complete examples are in ``examples/simple_api.py``, ``examples/classify_fast.py``,
and ``examples/classify_fractal.py``.

Legacy script note
------------------

``examples/orbits_and_cuts.py`` uses ``pysimple.trace_orbit()`` and computes a
simple toroidal-plane cut in Python. Prefer ``examples/simple_api.py`` when you
just want the shortest supported entry point.

Testing
-------

See ``test/python`` for pytest-based validation of the API. Activate ``.venv``
first so the tests use the same interpreter that has ``pysimple`` and its
dependencies installed.
