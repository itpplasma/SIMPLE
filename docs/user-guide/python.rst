Python interface
================

``pysimple`` is a thin Python layer over SIMPLE's Fortran library interface
``src/simple_api.f90``. It runs the same code as ``simple.x``: a Python run
configured from a ``simple.in`` file reproduces ``simple.x`` exactly (this is
tested in ``test/python/test_simple_api.py``). Installation is described in
:doc:`/getting-started/installation`.

Typical use:

.. code-block:: python

   import pysimple

   pysimple.init(config="simple.in")                  # exactly as simple.x
   pysimple.init("wout.nc", trace_time=1e-3)          # file and settings
   pysimple.init("wout.nc", config="simple.in", n_d=2)  # file and overrides

   start = pysimple.sample_surface(256, s=0.6)        # (5, n) array
   result = pysimple.trace_parallel(start)            # OpenMP over particles
   pysimple.write_output()                            # simple.x output files

``init`` starts from the compiled defaults, then applies a config file, an
equilibrium file and keyword settings, in this order. Keyword names are
``&config`` names and are case-insensitive; unknown names raise
``ValueError``. ``pysimple.params`` exposes derived values such as
``params.v0``.

States are :math:`(s, \vartheta, \varphi, v/v_0, v_\parallel/v)` in reference
coordinates with shape ``(5, n)``. ``trace_parallel`` returns a dictionary
with ``final_positions``, ``loss_times``, ``trap_parameter`` and
``perpendicular_invariant``. ``loss_times`` is ``trace_time`` for confined,
``-1`` for skipped and ``NaN`` for numerically unresolved particles.

Further functions trace single orbits (``trace_orbit``), compute bounce and
precession frequencies (``compute_canonical_frequencies``), trace to the next
banana tip or field-period plane (``trace_to_cut``), and convert between
states and invariants (``invariants_from_state``,
``states_from_invariants``). The full overview is the package README,
:doc:`/overview`, and the reference is
:doc:`/api/index`. :doc:`/tutorials/optimisation-loop` uses ``pysimple``
inside an optimisation loop.

.. note::

   At this revision, ``pysimple.sample()`` and ``pysimple.sample_surface(n, s)``
   with ``s`` equal to the ``sbeg`` given to ``init`` can return start points
   from uninitialised memory: ``simple_sample`` reallocates the start-surface
   arrays (``src/simple_api.f90``, ``allocate_particles``) but skips the
   field-line integration that fills them when the surface has not changed
   (``use_sampling_field``). Check that ``start[0]`` equals the requested
   surface, as :doc:`/tutorials/optimisation-loop` does, or use a surface
   different from ``sbeg`` in ``sample_surface``.

Fortran programs use the same interface through ``use simple_api``: set
``params`` variables or call ``simple_read_config``, then ``simple_init``,
``simple_sample``, ``simple_trace`` and ``simple_get_results``.
