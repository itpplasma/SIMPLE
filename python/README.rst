pysimple
========

``pysimple`` is a thin Python layer over SIMPLE's Fortran library interface
(``src/simple_api.f90``). It runs the same code path as ``simple.x``; a Python
run configured from a ``simple.in`` file reproduces ``simple.x`` exactly.

Installation
------------

.. code-block:: bash

    ./setup-venv.sh
    source .venv/bin/activate
    python -c "import pysimple"

Configuration
-------------

A run is configured by the ``&config`` namelist of ``simple.in``. The complete
list of keys, with defaults and meaning, is ``examples/simple_full.in``.
``init`` starts from the defaults, then applies a config file, an equilibrium
file and keyword settings, in this order:

.. code-block:: python

    import pysimple

    pysimple.init(config="simple.in")                      # as simple.x
    pysimple.init("wout.nc", trace_time=1e-3, ntestpart=32)  # no file
    pysimple.init("wout.nc", config="simple.in", n_d=2)      # file + overrides

Keyword names are ``&config`` names and are case-insensitive. Unknown names
raise ``ValueError``. The equilibrium argument sets ``netcdffile``,
``field_input`` and ``coord_input``. A Boozer chartmap selects
``isw_field_type=2`` unless the field type is given. Read derived values from
``pysimple.params``, for example ``pysimple.params.v0``.

Tracing
-------

States are ``(s, theta, phi, v/v0, v_par/v)`` in reference coordinates, with
shape ``(5, n)``. Angles are radians and are not wrapped.

.. code-block:: python

    start = pysimple.sample()                  # as configured by startmode
    start = pysimple.sample_surface(32, s=0.3)
    start = pysimple.sample_volume(32, s_inner=0.1, s_outer=0.9)
    start = pysimple.load_particles("start.dat")

    result = pysimple.trace_parallel(start)    # integrator="midpoint" overrides integmode
    pysimple.write_output()                    # simple.x output files

``loss_times`` is ``trace_time`` for confined particles, ``-1`` for particles
skipped by ``contr_pp`` (deep passing; set ``contr_pp=-1e10`` to trace them)
and ``NaN`` for numerical failures. ``classify_parallel``, ``classify_fast``
and ``classify_fractal`` add the trapped/passing flag and the J_parallel,
topological and fractal classes.

Single orbits
-------------

``trace_orbit(state, return_trajectory=True)`` returns the ``ntimstep``
samples of one orbit.

``compute_canonical_frequencies(state, n_periods=1)`` returns the
bounce/transit period, the toroidal shift per period and the frequencies
``omega_b`` and ``omega_phi`` in rad/s. A trapped cycle runs between
same-direction ``v_par = 0`` crossings, a passing cycle between signed
``2*pi`` poloidal advances. Check ``status`` against ``FREQ_*``; a lost orbit
reports ``FREQ_ORBIT_LOST``. Species and energy follow ``&config``
(``n_e``, ``n_d``, ``facE_al``); the default is a 3.5 MeV alpha.

``trace_to_cut(state, cut="tip")`` traces to the next ``v_par = 0`` point,
``cut="toroidal"`` to the next field-period plane. It needs a symplectic
integrator.

``invariants_from_state(state)`` returns ``h0 = (v/v0)**2``,
``j_perp = h0 * (1 - lambda**2) / B`` and ``p_phi = (c/q) P_phi`` in
magnetic-flux units (zero on axis, outward positive in its flux term; not the
radial label ``s``). ``axis_flux_native`` and ``edge_flux_native`` are
``A_phi`` at the axis and at the LCFS.

``states_from_invariants(h0, j_perp, p_phi)`` inverts this for axisymmetric
fields on the section ``grad(psi) x grad(B) = 0``. It returns all regular
roots on both ``v_par`` branches, ordered by major radius from the high- to the
low-field side (POTATO's ``R_c``). ``states_from_potato_invariants`` first
converts POTATO's ``psi_star`` with the given ``psi_axis``, ``psi_edge`` and
``v0_ratio = v0_POTATO / v0_SIMPLE``.

Fortran
-------

Fortran programs use the same interface through ``use simple_api``: set
``params`` variables or call ``simple_read_config``, then ``simple_init``,
``simple_sample``, ``simple_trace`` and ``simple_get_results``. ``simple.x``
itself is the sequence in ``simple_main::main``.
