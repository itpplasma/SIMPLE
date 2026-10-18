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

IMAS IDS workflow
-----------------

Install ``pysimple[imas]`` for the optional IMAS-Python interface and libneo's
equilibrium converter. ``pysimple.imas`` reads real Data Dictionary 3 and 4
objects; the default output dictionary is 4.1.1. A complete frozen-equilibrium
cohort workflow is available in memory or through an IMAS DBEntry::

    from pysimple.imas import run_entry

    result = run_entry(
        "input.nc", "output.nc", duration=1e-3, source_interval=0.01,
        core_profiles=True, wall=True,
    )

The command-line equivalent is::

    simple-imas input.nc output.nc --duration 1e-3 --source-interval 0.01 \
        --core-profiles --wall

Plain ``.nc`` paths use IMAS-Python's portable NetCDF backend. Standard DBEntry
URIs can select other available backends. Input contains ``equilibrium`` and
``distribution_sources``; ``--marker-ids distributions`` selects population
markers instead. ``core_profiles`` and ``wall`` are optional. Existing output
entries are not overwritten. ``--time`` selects an exact common initial slice;
the field and bulk plasma remain frozen throughout the run.

``run(markers, equilibrium=eq, duration=..., ...)`` accepts IDS objects directly.
An existing VMEC or Boozer chartmap can replace ``equilibrium`` with ``field``
for collisionless marker workflows. The native configuration remains available
through keyword arguments or a JSON ``--settings`` file. Configuration that
sets species, physical input, normalization or duration belongs to the actor.

Input and output conventions:

* Guiding-centre ion markers use coordinate identifiers 4, 5, 3, 301 and 403:
  R [m], phi [rad], Z [m], kinetic energy [eV], and v_parallel/v. Column order
  may vary. Multiple ion populations are traced separately with their own
  mass number and charge. SIMPLE uses its native mass-number times proton-mass
  convention. Actual-particle markers and electron populations are rejected.
* Axisymmetric equilibria require a uniform rectangular R/Z grid, poloidal
  flux, F, pressure and q profiles, the magnetic axis and a closed boundary.
  DD3 COCOS11 and DD4 COCOS17 flux signs are mapped explicitly to libneo's
  COCOS3 intermediary. The converter runs in an isolated child process.
  Conversion resolution must be checked for the equilibrium being studied.
* Collision profiles use electron temperature and density plus one or two
  bulk-ion densities and temperatures in eV and m^-3. IDS rho_tor_norm is
  squared for the native collision abscissa s. Quasineutrality is checked;
  rotation, mixed charge states and extra species are rejected. The current
  native collision operator supports alpha markers (A=4, charge=2) only.
* Walls accept Cartesian triangle GGD meshes with standard coordinate
  identifiers and one-based node indices, or one closed R/Z limiter outline
  revolved into an STL. Lengths are metres. Tracing walls requires a build
  with ``SIMPLE_ENABLE_CGAL=ON``. Toroidal discretization of a revolved wall
  needs convergence checks for wall-load studies.
* ``distribution_sources`` marker weights are particles/s. An explicit
  ``source_interval`` converts these rates to the number of particles in a
  finite cohort. ``distributions`` input already has particle weights and
  rejects a source interval. These conventions are never interchanged silently.

Output ``distributions`` contains initial and final marker snapshots, particle
counts, and total kinetic energy [J] for every population. Lost markers have
zero final weight. Private negative coordinate identifiers -1 and -2 retain
termination time and native exit condition (0 completed, 1 LCFS, 2 wall).
Numerical failures and skipped markers raise an error. Snapshots describe the
finite cohort; they are not a continuous-source residence-time histogram.

ASCOT5 compatibility follows its current
``a5py/templates/imasinterface.py`` marker identifiers, SI/eV/radian conventions
and Cartesian GGD node indexing. ASCOT's HDF5 marker API uses degrees internally;
its IDS reader performs that conversion. SIMPLE keeps radians. ASCOT and SIMPLE
have different solver models and output-distribution semantics; compare matched
physics, actual end times, boundaries and weights before comparing results.

Fortran
-------

Fortran programs use the same interface through ``use simple_api``: set
``params`` variables or call ``simple_read_config``, then ``simple_init``,
``simple_sample``, ``simple_trace`` and ``simple_get_results``. ``simple.x``
itself is the sequence in ``simple_main::main``.
