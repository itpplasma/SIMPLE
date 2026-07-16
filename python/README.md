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

    n_lost = ((results["loss_times"] > 0.0) & (results["loss_times"] < 5e-5)).sum()
    n_skipped = (results["loss_times"] < 0.0).sum()
    print(f"Lost particles: {n_lost} (skipped deep-passing: {n_skipped})")

Complete examples are in ``examples/simple_api.py``, ``examples/classify_fast.py``,
and ``examples/classify_fractal.py``.

Canonical orbit frequencies
---------------------------

``compute_canonical_frequencies`` traces one particle until it completes the
requested number of poloidal cycles:

.. code-block:: python

    position = [0.6, 0.0, 0.25 * np.pi, 1.0, 0.3]
    frequency = pysimple.compute_canonical_frequencies(
        position, integrator="midpoint", n_periods=4
    )
    print(frequency["omega_b"], frequency["omega_phi"])

For trapped particles, a cycle is measured between same-direction
``v_parallel = 0`` crossings. For passing particles, it is measured between
successive signed ``2*pi`` advances of the unwrapped poloidal angle. The
toroidal displacement is also unwrapped. Periods are seconds, displacements
are radians, and both frequencies are angular frequencies in rad/s.

Use ``n_periods=1`` for the inexpensive single-cycle value in an axisymmetric
field. Larger values return the mean and sample spread over several cycles,
which is useful in asymmetric fields. Always inspect ``status`` before using a
result; losses, integrator errors, and the maximum-step limit remain visible
instead of being silently discarded. An orbit that leaves the plasma through
the outer radial boundary reports ``FREQ_ORBIT_LOST``.

The particle species and energy follow the same defaults as ``simple.in``
(3.5 MeV alphas). A frequency computed for a different particle needs the
matching overrides in ``init``, for example a 5 keV deuteron:

.. code-block:: python

    pysimple.init("field.nc", deterministic=True, n_e=1, n_d=2, facE_al=700.0)

A fast particle in a small equilibrium can be genuinely unconfined; a
``FREQ_ORBIT_LOST`` result for a seed that should be confined usually means
the species/energy overrides are missing.

The same event machinery is available as ``pysimple.trace_to_cut``. Select
``cut="tip"`` for a same-direction ``v_parallel=0`` section or
``cut="toroidal"`` for a field-period section. The returned state is converted
back to public reference coordinates. This driver currently requires one of
SIMPLE's symplectic integrators.

Legacy script note
------------------

``examples/orbits_and_cuts.py`` uses ``pysimple.trace_orbit()`` and computes a
simple toroidal-plane cut in Python. Prefer ``examples/simple_api.py`` when you
just want the shortest supported entry point. The single-particle state is
``[s, theta, phi, v/v0, lambda]`` with ``lambda = v_parallel / v``; the example
includes both a confined trapped seed and a confined passing seed because the
sign of ``lambda`` alone does not determine the orbit topology. The script also
sets ``contr_pp=-1e10`` so deep-passing seeds are traced instead of being
skipped by the default ``contr_pp=-1`` policy.

Testing
-------

See ``test/python`` for pytest-based validation of the API. Activate ``.venv``
first so the tests use the same interpreter that has ``pysimple`` and its
dependencies installed.
