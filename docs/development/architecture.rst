Architecture
============

This page maps the code: how a run flows through the modules, where fields,
coordinates and integrators live, and what to change to add a field source or
an integrator. Implementation detail of coordinates and fields is in
:doc:`design/coordinates-and-fields`; the equations are in :doc:`/model/index`.

Run sequence
------------

``simple.x`` calls ``simple_main::main`` (``src/simple_main.f90``):

1. ``load_config`` reads ``&config`` into module ``params`` and ``&profiles``
   into ``simple_profiles``.
2. ``init_run`` loads the field (``init_field``), derives the normalisation
   (``params_init``: :math:`v_0`, :math:`\rho_0`, step :math:`\Delta\tau`) and
   initialises collisions.
3. ``init_sampling_field`` and ``sample_start`` fill ``zstart`` according to
   ``startmode`` (``src/samplers.f90``).
4. ``init_tracing_field`` selects the integration field.
5. ``trace_parallel`` runs ``trace_orbit`` for each marker in an OpenMP loop
   (``schedule(dynamic, 1)``); with ``gpu_mode = 'production'``,
   ``trace_gpu_production`` replaces it.
6. ``write_output`` writes the text files and ``results.nc``.

``simple_api`` (``src/simple_api.f90``) exposes the same steps as separate
calls (``simple_init``, ``simple_sample``, ``simple_trace``, ...) plus
single-orbit diagnostics. Public states are
:math:`(s, \vartheta, \varphi, v/v_0, \lambda)` in reference coordinates;
``trace_orbit`` converts them with ``ref_to_integ`` before integration.

Module map
----------

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Module (file)
     - Role
   * - ``params`` (``src/params.f90``)
     - ``&config`` keys, run state, exit codes, per-marker random streams
   * - ``simple_main``, ``simple_api``
     - run sequence, per-orbit loop, wall and LCFS events, recovery, output
   * - ``simple`` (``src/simple.f90``)
     - ``tracer_t``, normalisation, ``encode_symplectic_state``
   * - ``samplers``, ``fusion_source``
     - start positions; fusion-weighted births
   * - ``field`` (``src/field.F90``)
     - ``field_from_file``: field object from ``field_input``
   * - ``src/field/field_*.f90``
     - field classes (coils, Boozer chartmap, splined) and canonical
       fields (``field_can_flux``, ``_boozer``, ``_meiss``, ``_albert``,
       ``_spectre``, ``_test``)
   * - ``field_can_mod`` (``src/field_can.f90``)
     - canonical field dispatch, ``integ_to_ref`` / ``ref_to_integ``
   * - ``magfie_sub`` (``src/magfie.f90``)
     - non-canonical field evaluation for RK45 and sampling
   * - ``reference_coordinates``, ``src/coordinates/``
     - reference coordinate system ``ref_coords``, radial scaling
   * - ``orbit_symplectic_base``, ``orbit_symplectic``,
       ``orbit_symplectic_quasi``, ``orbit_symplectic_euler1``
     - symplectic and quasi-symplectic steppers, LCFS event location
   * - ``alpha_lifetime_sub`` (``src/sub_alpha_lifetime_can.f90``)
     - RK45 guiding-centre equations (``velo_can``) and recovery steps
   * - ``orbit_fo_boris``, ``orbit_fo_field``
     - full-orbit Boris pusher
   * - ``spectre_orbit``, ``spectre_sympl_orbit``, ``spectre_fo_field``,
       ``spectre_fo_hybrid``, ``field_line_spectre``
     - SPECTRE tracing per volume, full orbit near interfaces, field lines
   * - ``collis_alp`` (``src/collis_alphas.f90``)
     - Coulomb collision operator (``stost``)
   * - ``classification``, ``detect_oneline_mod``, ``cut_detector``
     - orbit classification
   * - ``stl_wall_intersection``, ``src/wall/stl_wall_cgal.cpp``
     - STL wall segment test through CGAL
   * - ``orbit_frequencies``, ``orbit_invariants``
     - bounce and toroidal frequencies; invariants to starting states
   * - ``netcdf_results_output``, ``netcdf_orbit_output``
     - ``results.nc``, ``orbits.nc``
   * - ``simple_gpu``, ``simple_cuda_native``, ``*.cu``
     - OpenACC and native CUDA tracing for Boozer fields
   * - ``python/pysimple/``
     - Python layer over ``simple_api`` and ``params``

Libraries: libneo supplies VMEC splines, Boozer conversion, coordinate
systems (``libneo_coordinates``), the field base class ``magnetic_field_t``,
Biot–Savart, SPECTRE fields and the ODE and spline routines; fortplot draws
diagnostic plots; fortnum provides the root finder and the GPU RK5(4)
integrator. Versions are pinned in
``cmake/dependency_versions.json``.

Fields and coordinates
----------------------

Three layers separate field data from the coordinates used for integration.

1. **Field objects** extend libneo ``magnetic_field_t`` and implement
   ``evaluate(x, Acov, hcov, Bmod)`` in their own coordinates: ``vmec_field_t``
   (libneo), ``coils_field_t``, ``boozer_chartmap_field_t``,
   ``splined_field_t`` (any field resampled on reference coordinates) and
   ``spectre_field_t`` (libneo). ``field_from_file`` chooses the class from
   the file.
2. **Reference coordinates** (``ref_coords``, a libneo
   ``coordinate_system_t``) are VMEC :math:`(s,\vartheta,\varphi)` or a
   chartmap :math:`(\rho,\vartheta,\zeta)`, set by ``coord_input``.
3. **Integration coordinates** are chosen by ``integ_coords``
   (:doc:`/user-guide/field-sources`). Canonical choices fill a
   ``field_can_t`` with :math:`A_\vartheta, A_\varphi, h_\vartheta,
   h_\varphi, B` and their first and second derivatives, then add
   :math:`H, p_\vartheta, v_\parallel` and their derivatives in
   :math:`(r,\vartheta,\varphi,p_\varphi)`. ``field_can_from_name`` sets the
   procedure pointers ``evaluate``, ``integ_to_ref`` and ``ref_to_integ``.
   RK45 uses ``magfie``, a procedure pointer set by ``init_magfie(id)``, which
   returns :math:`B`, :math:`\sqrt g`, :math:`\nabla\ln B`, :math:`h_i`,
   :math:`h^i` and :math:`(\nabla\times\mathbf h)^i`.

Integrators
-----------

``tracer_t`` (``src/simple.f90``) holds one marker's ``field_can_t``,
``symplectic_integrator_t``, ``multistage_integrator_t`` (Gauss and Lobatto
stages) and full-orbit state. ``orbit_sympl_init`` maps ``integmode``
(:doc:`/user-guide/running`) to a step routine ``raw_timestep_sympl`` and its
quasi-symplectic fallback; ``orbit_timestep_sympl`` wraps it with the LCFS
event search. RK45 (``integmode = 0``) integrates ``velo_can`` with libneo
``odeint_allroutines``. The GPU paths do not use procedure pointers: they
evaluate the Boozer field locally and reuse the step algebra of
``orbit_symplectic_euler1`` and ``orbit_rk_core``.

Python binding
--------------

``python/CMakeLists.txt`` wraps ``simple_api.f90``, ``params.f90``,
``orbit_symplectic_base.f90``, ``orbit_frequencies.f90`` and
``orbit_invariants.f90`` with f90wrap into the extension
``_simple_backend``. ``pysimple`` adds array handling and plotting only;
physics stays in Fortran so that ``simple.x`` and ``pysimple`` agree. A
Fortran program can link the CMake target ``simple`` and ``use simple_api``.

Adding a field source
---------------------

1. Write ``type, extends(magnetic_field_t)`` in ``src/field/`` with
   ``evaluate`` returning covariant :math:`A_i`, :math:`h_i` and :math:`B` in
   Gaussian units (G, cm) in its native coordinates; convert units on load as
   ``create_coils_field`` does.
2. Detect the file in ``field_from_file`` (``src/field.F90``). If the native
   coordinates are not the reference coordinates, wrap the field with
   ``create_splined_field`` as for coils; it then runs with Meiss, Albert or
   RK45 on ``REFCOORDS``.
3. A field that needs its own canonical coordinates also needs a
   ``field_can_<name>.f90`` with ``evaluate``, ``integ_to_ref`` and
   ``ref_to_integ``, a case in ``field_can_from_name`` and an id in
   ``magfie_sub``.
4. Add a roundtrip test and an orbit comparison against an existing source on
   the same equilibrium, and update :doc:`design/coordinates-and-fields`.

Adding an integrator
--------------------

1. Add a constant in ``orbit_symplectic_base`` and a step routine in
   ``orbit_symplectic.f90`` that advances ``si%z`` from a ``field_can_t``
   and returns a ``SYMPLECTIC_STEP_*`` status; reuse the Newton and LCFS
   event helpers.
2. Add the case to ``orbit_sympl_init`` and, if it has stages, size
   ``multistage_integrator_t``.
3. Document the new ``integmode`` value in ``examples/simple_full.in``.
4. Test the convergence order on the analytic field (``integ_coords = -1``),
   the boundedness of :math:`H` and :math:`p_\varphi` and a short loss
   comparison against an existing scheme.

Invariants to keep
------------------

- Gaussian CGS inside the code; SI and eV only at input and output.
- Angles in radians, never wrapped during integration.
- The public state is :math:`(s, \vartheta, \varphi, v/v_0, \lambda)` in
  reference coordinates; :math:`s` stays normalised toroidal flux at every
  interface, also for chartmaps where the coordinate is :math:`\rho`.
- Fields and canonical maps are built before the parallel region and read
  only inside it; per-marker state lives in ``tracer_t`` or ``threadprivate``
  variables.
- No procedure pointers or ``class()`` dispatch inside GPU kernels; the CPU
  path is unchanged when acceleration is off.

Subsystem invariants for tests are listed in :doc:`testing`. Known design
limitations are tracked in issues
`#288 <https://github.com/itpplasma/SIMPLE/issues/288>`_ (VMEC inversion),
`#290 <https://github.com/itpplasma/SIMPLE/issues/290>`_ (``sqgBctr``),
`#291 <https://github.com/itpplasma/SIMPLE/issues/291>`_ (radial scaling),
`#292 <https://github.com/itpplasma/SIMPLE/issues/292>`_ (global ``magfie``)
and `#293 <https://github.com/itpplasma/SIMPLE/issues/293>`_ (module-level
canonical coordinates).
