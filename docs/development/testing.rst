Testing
=======

Tests are registered with CTest and run through ``make`` from the repository
root. Every change of behaviour comes with a test that fails without it.

Test targets
------------

All targets except ``test-all`` first rebuild with deterministic floating
point (``SIMPLE_DETERMINISTIC_FP=ON``).

=========================  ==================================================
Command                    Runs
=========================  ==================================================
``make test-fast``         all but labels ``slow``, ``regression``,
                           ``performance``; run before every push
``make test``              all but ``regression``
``make test-nopy``         as ``make test`` without label ``python``
``make test-slow``         label ``slow``
``make test-smoke``        label ``smoke`` (sub-minute gate used by libneo)
``make test-regression``   label ``regression``: golden records vs ``main``
``make test-all``          every test, current build configuration
=========================  ==================================================

``TEST=<regex>`` selects tests by name, ``VERBOSE=1`` shows all output and
``CONFIG=Fast`` selects the build profile. The regex passes through ``make``:
quote it and write ``$`` as ``$$``, e.g. ``TEST='^test_finite_check$$'``.
Tests labelled ``python`` and some unlabelled wrappers need the packages from
``./setup-venv.sh``.

Memory allocation on workstations and CI
---------------------------------------

Configure on the machine that will run the tests. CMake writes
``build/ctest-resources.json`` with one ``mem_gb`` resource whose slots each
represent one GiB. Usable memory is the minimum of Linux ``MemAvailable``,
cgroup memory headroom, and any positive Slurm memory allocation
(``SLURM_MEM_PER_NODE`` or ``SLURM_MEM_PER_CPU * SLURM_CPUS_ON_NODE``, in MiB).
For cgroup v2, CMake locates the job via ``/proc/self/cgroup``, takes the
smallest finite ``memory.max`` from that cgroup up to the mount root, and
subtracts the job cgroup's ``memory.current``. The v1 fallback uses
``memory.limit_in_bytes`` and ``memory.usage_in_bytes`` under the memory
controller, ignoring huge unlimited sentinels. Negative headroom becomes zero.
Missing cgroup files impose no constraint; missing ``MemAvailable`` falls back
to CMake's available physical memory query.
The automatic budget is this usable memory, rounded down to GiB, less a
reserve of 25%, rounded up, or 4 GiB, whichever is larger. This leaves room
for the OS, other jobs and tests without reservations. Cap this snapshot on
a shared machine, and reconfigure when its available memory changes::

    cmake -S . -B build -DSIMPLE_TEST_MEMORY_GB=24
    ctest --test-dir build/test -j8 --output-on-failure

CMake records ``CTEST_RESOURCE_SPEC_FILE`` in every generated
``CTestTestfile.cmake``. Both ``ctest`` from ``build`` and
``ctest --test-dir build/test`` therefore load the generated spec automatically.
Thus CI's existing ``ctest -j16`` also uses memory allocation without another
command-line option. A site-specific spec can override the default::

    ctest --test-dir build/test -j16 --resource-spec-file /path/to/resources.json

For example, a 24 GiB budget has this format::

    {"version": {"major": 1, "minor": 0},
     "local": [{"mem_gb": [{"id": "ram", "slots": 24}]}]}

CTest's ``RESOURCE_GROUPS`` alone does **not** limit parallelism without a
resource spec; tests still run. This is documented in the
`CTest resource allocation manual
<https://cmake.org/cmake/help/latest/manual/ctest.1.html#resource-allocation>`_
and verified with CTest 4.4.4. The automatic default above prevents that
unbounded fallback. A test requiring more slots than the budget is reported
as ``Not Run`` with insufficient resources; the budget is never raised to
force it to run. Memory reservations govern scheduling within one CTest
invocation, not OS memory limits or separate CTest invocations. Keep OpenMP
thread counts appropriate to the CPU allocation as well.

Reservations are peak RSS rounded up to whole GiB, measured one registered
test at a time with GNU ``/usr/bin/time -v`` inside
``systemd-run --user --scope -p MemoryMax=16G``. The workstation measurements
use GNU Fortran 16.2.1, the deterministic ``Fast`` build, four OpenMP threads
and one OpenBLAS thread; golden-record scripts retain their fixed two threads
and the wall test retains its fixed single thread.

==============================================  ==========  =====
Test                                            RSS (GiB)   Slots
==============================================  ==========  =====
test_splined_field_derivatives                  0.8202      1
test_bder_from_splines                          0.8191      1
test_magfie_refcoords_splines                   1.1990      2
test_hcurl_from_dhcov_splines                   0.8194      1
test_splined_kernel_bench                       12.6803     13
test_chartmap_pipeline                          9.8633      10
test_chartmap_rz_consistency                    9.8601      10
gc_ck_ncsx_vmec_vs_meiss                        11.9712     12
golden_record_albert_coils                      14.0856     15
golden_record_boozer_ncsx                       7.2261      8
golden_record_meiss_coils                       12.0610     13
test_collisional_stl_wall                       10.3197     11
generate_chartmap_ncsx_map2disc                 0.3585      1
generate_chartmap_ncsx_map2disc_wall            0.3940      1
test_chartmap_ncsx_map2disc_wall_offset         0.3940      1
test_chartmap_boundary_param                    0.3988      1
==============================================  ==========  =====

Re-measure when adding NCSX consumers or increasing fixture sizes.
Golden cases use ``flock`` (util-linux) only while preparing their shared
reference build and input files; their numerical runs share the memory budget.

Kinds of tests and labels
-------------------------

=========================================  =====================================
Kind (labels)                              Location: checks
=========================================  =====================================
Unit (``unit``)                            ``test/tests/test_*.f90``: one
                                           procedure against an independent
                                           oracle
Integration (``integration``, ``system``)  ``test/tests/``: several modules or
                                           ``simple.x`` on a small case;
                                           physical invariants
Python (``python``)                        ``test/python/test_*.py`` (pytest):
                                           ``pysimple``, ``tools/``, agreement
                                           of ``pysimple`` and ``simple.x``
Golden record (``regression``)             ``test/golden_record/<case>/``:
                                           end-to-end output against ``main``
=========================================  =====================================

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

=======================  ===================================================
Subsystem                Must hold after a change
=======================  ===================================================
Coordinate transforms    ``ref -> integ -> ref`` roundtrip, periodicity,
                         explicit :math:`s`, :math:`\sqrt{s}` or :math:`\rho`
                         at every interface
VMEC, Boozer, chartmap   same start surface, same short-time orbit within
                         tolerance, statistically equal loss fraction
Symplectic integrators   bounded invariants, no drift beyond the reference
                         envelope
Collisions               Maxwellian fixed point, analytic pitch and energy
                         relaxation, seeded reproducibility
Wall losses              same wall hit in Cartesian and coordinate space
NetCDF output            unchanged meaning, units and frame of ``zend``,
                         ``xend_cart``, ``wall_hit_cart``
GPU and OpenMP           CPU path unchanged when acceleration is off,
                         short-run CPU/GPU agreement, long-run statistical
                         agreement, no CPU slowdown
=======================  ===================================================

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
