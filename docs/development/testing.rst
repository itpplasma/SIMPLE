Testing
=======

Tests are registered with CTest and run through ``make`` from the repository
root. Every change of behaviour comes with a test that fails without it.

Test targets
------------

All targets except ``test-all`` first rebuild with deterministic floating
point (``SIMPLE_DETERMINISTIC_FP=ON``).

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Command
     - Runs
   * - ``make test-fast``
     - everything except labels ``slow``, ``regression``, ``performance``;
       run before every push
   * - ``make test``
     - everything except ``regression``
   * - ``make test-nopy``
     - as ``make test`` without label ``python``
   * - ``make test-slow``
     - label ``slow`` only
   * - ``make test-smoke``
     - label ``smoke`` (sub-minute gate used by libneo)
   * - ``make test-regression``
     - label ``regression``: golden records against ``main``
   * - ``make test-all``
     - every test, in the current build configuration

``TEST=<regex>`` selects tests by name, ``VERBOSE=1`` shows all output and
``CONFIG=Fast`` selects the build profile. The regex passes through ``make``:
quote it and write ``$`` as ``$$``, e.g. ``TEST='^test_finite_check$$'``.
Tests labelled ``python`` and some unlabelled wrappers need the packages from
``./setup-venv.sh``.

Kinds of tests and labels
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 33 45

   * - Kind (label)
     - Location
     - Checks
   * - Unit (``unit``)
     - ``test/tests/test_*.f90``
     - one procedure against an independent oracle
   * - Integration, system (``integration``, ``system``)
     - ``test/tests/``
     - several modules or ``simple.x`` on a small case; physical invariants
   * - Python (``python``)
     - ``test/python/test_*.py`` (pytest)
     - ``pysimple``, ``tools/``, agreement of ``pysimple`` and ``simple.x``
   * - Golden record (``golden_record``, ``regression``)
     - ``test/golden_record/<case>/simple.in``
     - end-to-end output against ``main``

Add ``slow`` to tests that run longer than about a minute. Register Fortran
tests in ``test/tests/CMakeLists.txt`` and pytest files with ``add_pytest`` in
``test/python/CMakeLists.txt``.

Oracles and tolerances
----------------------

- The oracle is independent of the implementation: an analytic result, an
  identity, a roundtrip (``ref -> integ -> ref``), a symmetry, a convergence
  rate or a second field representation. Rerunning the same code is not an
  oracle.
- A tolerance follows from a physical or numerical scale (``relerr``, spline
  order and grid spacing, Monte Carlo error :math:`\sim 1/\sqrt{N}`), never
  from the observed error.
- Stochastic results are compared as distributions or moments, not particle
  by particle.
- Data available on only one machine must not be required; missing optional
  data skips the test.

Invariants to keep
------------------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Subsystem
     - Must hold after a change
   * - Coordinate transforms
     - ``ref -> integ -> ref`` roundtrip, periodicity, explicit
       :math:`s`, :math:`\sqrt{s}` or :math:`\rho` at every interface
   * - VMEC, Boozer, chartmap
     - same start surface, same short-time orbit within tolerance,
       statistically equal loss fraction
   * - Symplectic integrators
     - bounded invariants, no drift beyond the reference envelope
   * - Collisions
     - Maxwellian fixed point, analytic pitch and energy relaxation, seeded
       reproducibility
   * - Wall losses
     - same wall hit in Cartesian and coordinate space
   * - NetCDF output
     - unchanged meaning, units and frame of ``zend``, ``xend_cart``,
       ``wall_hit_cart``
   * - GPU and OpenMP
     - CPU path unchanged when acceleration is off, short-run CPU/GPU
       agreement, long-run statistical agreement, no CPU slowdown

Golden records
--------------

``golden_record_<case>`` clones and builds ``main`` under ``build/``, runs
``test/golden_record/<case>/simple.in`` with both builds and compares the
outputs with relative tolerance ``1.5e-7``. The first run needs network
access and several minutes. Use ``make test-regression TEST=golden_record_``:
it pins the reference libneo revision (``libneo_golden`` in
``cmake/dependency_versions.json``) for both builds. Patches in
``test/golden_record/reference_patches/`` make old reference revisions
buildable and comparable. The stochastic ``collision`` case uses its own
acceptance rule (:doc:`/collision-golden-contract`). Add a case when a new
field type, integrator or output path becomes part of the supported workflow.
A change of golden output is a scientific change and must be explained in the
pull request (:doc:`/contributing`).

The bounded numerical verification report and its evidence table are in
:doc:`/verification`.
