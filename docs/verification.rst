Numerical verification
======================

``tools/generate_verification_report.py`` runs a bounded subset of SIMPLE's
verification evidence and writes a Markdown report with PNG figures, raw data
and logs. It contributes evidence for the requested EUROfusion Standard
Software v2.2 criteria 3.1 and 3.4. It does not establish compliance with the
whole standard or validate the guiding-center model against experiment.

Existing evidence
-----------------

The source files below define the numerical contracts. Inclusion in this
inventory does not mean that the generated report runs every test.

.. list-table:: Evidence inventory
   :header-rows: 1
   :widths: 27 51 22

   * - Source
     - Reference and acceptance rule
     - Bounded report
   * - ``test/golden_record/``
     - Revision comparisons of Boozer, canonical, Meiss, Albert, collision,
       tokamak and classifier runs. ``compare_files.py`` defaults to
       rtol=1e-7, atol=1e-12 (CTest sets rtol=1.5e-7); classifier final states
       are excluded from that comparator. ``compare_orbit_results.py`` also
       checks physical validity of accepted recovery transitions. References
       are built with the pinned golden libneo revision. These are regression
       comparisons, not analytic oracles.
     - Not run: requires a second build.
   * - ``test/golden_record_sanity/``
     - Matching, mismatching and NaN fixtures test the comparison machinery.
       They do not verify an orbit calculation.
     - Not run.
   * - ``tools/step_convergence.py`` and
       ``test/python/test_step_convergence.py``
     - Paired outcomes at two steps: bias <= binomial SE in magnitude and
       flip rate <= 3 SE. KS and Wasserstein distances are descriptive;
       they are not part of the implemented verdict. Unit tests use synthetic
       ``times_lost.dat`` files, including survivors and known flips.
     - Not run: no physical order follows from those unit tests.
   * - ``test/tests/test_sympl.f90`` (``test_sympl_tokamak``)
     - Analytic circular tokamak; Hamiltonian spread over about 30 bounces.
       Euler1/Euler2/midpoint/Gauss4 relative bounds: 0.124/0.120/0.0053/4.5e-5.
     - New short-time midpoint study uses the same analytic field backend.
   * - ``test/tests/test_axis_crossing.f90``
     - Six near-axis starts at two resolutions; energy excursion <= 1e-3/1e-2
       and half-trace mean-energy difference <= 1e-6/1e-4, with no axis loss.
     - Not run: longer integration.
   * - ``test/tests/test_orbit_invariants.f90``
     - Analytic forward definitions, inverse reconstruction of H0/J_perp/p_phi,
       canonical normalization and POTATO flux conventions. Reconstruction
       scaled tolerances: 1e-11/1e-10/1e-9.
     - Run; unprinted residuals are reported as bounds only.
   * - ``test/tests/test_profiles.f90``
     - Flat/scalar collision coefficients agree to relative 1e-12; seeded
       thermalized fractions and final coordinates agree. Slowing-down
       residence histogram matches p²/(p³+pc³) within 4%; Spitzer drag within
       3% with 0.2 < pc/v0 < 0.5. Maxwellian fixed point: 8000 markers, 4000
       collision steps, final/equilibrium mean p² within 5%, populated-bin
       deviation from initial and analytic histograms below 10%.
     - Run; printed distribution and moment errors are retained.
   * - ``test/tests/test_boozer_chartmap.f90`` and
       ``test_chartmap_scaling.f90``
     - Analytic fixture: A_theta=torflux*rho², zero radial h, rho/s roundtrip
       <= 1e-14. Field/length/flux scaling, interpolation and Cartesian inverse
       checks use their existing tolerances.
     - Run; assertions bound unprinted residuals.
   * - ``test/tests/test_boozer_chartmap_roundtrip.f90``
     - VMEC to Boozer export/reimport; 50 field points with relative B error
       <= 8e-5; 500-step direct/file midpoint orbits with absolute s error
       <= 2e-9. Both paths share the equilibrium and implementation.
     - Run; both measured errors are retained.
   * - ``test/tests/test_chartmap_rz_consistency.f90``
     - R/Z inverse and VMEC/chartmap geometry comparisons, including boundary
       agreement within 1.5 cm.
     - Not run.
   * - ``test/tests/test_spectre_sympl_volume.py``
     - Fixture-based midpoint/Euler energy excursion and fitted secular
       growth below 5e-6; cross-scheme and short-time RK45 comparisons.
       Separate ``field_can/test_magfie_spectre_analytic.f90`` compares drift
       derivatives with central differences to relative 1e-6.
     - Not run.
   * - ``test/tests/test_fo_boris.f90`` and
       ``field_can/test_boozer_spline_derivatives.f90``
     - Full-orbit energy/adiabatic-invariant and geometry checks; spline first
       and second derivatives compared with finite differences (2e-6 and
       5e-3 relative tolerances, denominator floor 1e-10).
     - Not run.
   * - ``test/python/test_simple_api.py`` and ``test_orbit_macro.py``
     - API/executable parity, invariant reconstruction and matching sampled
       trajectories. Shared implementation parity is not independent physics.
     - Not run.

Regenerate and inspect
----------------------

Run from the repository root, with the documented compiler and library
prerequisites in ``README.md``. Python needs NumPy, Matplotlib and netCDF4;
pytest and SciPy are used by the report helper tests. The Python/Fortran
``pysimple`` bindings are not required for the report.

.. code-block:: console

   make build-deterministic CONFIG=Fast FLAGS='-DENABLE_PYTHON_INTERFACE=OFF'
   make test CONFIG=Fast TEST='^test_verification_report$$'
   python tools/generate_verification_report.py

Read ``build/output/verification-report/verification-report.md`` alongside its
local PNG images. The report directory also contains ``measurements.json``,
CSV measurements, ``consumer-tests.log``, and per-run namelists, starts,
NetCDF trajectories, exit codes and logs. Build outputs are ignored by Git.
Keep the directory together when sharing the Markdown.

The default subprocess budget is 300 seconds, including the selected
``make test`` call, so build first. Each invocation requires an empty output
directory inside ``build/`` to prevent stale evidence. To retain several runs:

.. code-block:: console

   python tools/generate_verification_report.py --output-dir build/output/verification-second

Missing tests, execution failures, nonfinite data, unresolved/skipped markers
and timeouts produce a nonzero exit code and an incomplete/failing report.
The selected existing checks run via ``make test``; experimental trajectories
run ``build/simple.x`` in isolated report subdirectories.

The optional ``verification_report=true`` input on the CI workflow dispatch
builds SIMPLE, tests the helper, generates the report and uploads its directory
even on failure. No report generation is added to ordinary push/PR CI; only
the small helper regression test is added to the existing test suite.

Interpreting the figures
------------------------

The convergence plot measures SIMPLE's midpoint endpoint error relative to a
256-step two-stage Gauss reference (mode 5, order four) on the analytic TEST tokamak. A 128-step Gauss run
must differ by less than 1% of the smallest midpoint error. The error is the maximum absolute difference in the two
unwrapped endpoint angles (radians) at normalized time 0.1. Gauss stores its
last-stage radius (``src/orbit_symplectic.f90``), with endpoint angular
updates; its radial and pitch outputs therefore cannot serve as endpoint
reference values in this study. Four steps (8, 16, 32, 64) must give fitted order between
1.8 and 2.2; all points enter the log-log fit. Relative energy excursion must
stay below 1e-3 on those short trajectories. Analytic definitions also measure
relative J_perp and absolute native toroidal canonical momentum excursions,
both bounded by 1e-10 in normalized TEST units (B0=R0=ro0=1, a=0.5).
These fixed diagnostic thresholds
are not production accuracy guarantees. RK45 is deliberately excluded from
this TEST-field comparison: ``src/magfie.f90`` documents its missing curvature
drift in that backend.

Sampling uncertainty uses three **measured** SIMPLE ensembles on the generated
analytic chartmap fixture, with independent pseudorandom starting angles and
pitch on rho=0.85. They describe that explicit sampling law and a 1 microsecond
trace, not a physical reactor loss prediction. Each seed supplies nested
prefixes of 32, 64, 128 and 256 markers; prefixes are correlated. The report
prints loss counts, binomial standard errors, 95% Wilson intervals and the
sample standard deviation/range across seeds at N=256. Skipped or unresolved
orbits invalidate a binomial estimate rather than being silently counted.

The binomial standard error sqrt(p(1-p)/N) and its worst-case 1/(2 sqrt(N))
are **calculated** sampling guidance. The `NIST Wilson score description
<https://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm>`_ specifies the
interval formula; helper tests compare against SciPy independently, including
zero and all losses. Intervals assume independent, equally weighted markers.
They exclude numerical bias, weighted/correlated sampling, equilibrium and
model uncertainty. Three seeds do not measure confidence-interval coverage.
