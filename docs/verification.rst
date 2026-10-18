Numerical verification report
=============================

This page maps SIMPLE's numerical claims to the test that checks each one,
its oracle and its tolerance, and explains how to generate and read the
verification report of ``tools/generate_verification_report.py``.

Evidence
--------

The report runs a bounded subset of the tests below (last column).

.. list-table::
   :header-rows: 1
   :widths: 26 52 22

   * - Source
     - Oracle and acceptance rule
     - In the report
   * - ``test/tests/test_sympl.f90`` (``test_sympl_tokamak``)
     - Analytic circular tokamak; relative spread of :math:`\bar H` over
       about 30 bounces below 0.124 / 0.120 / 0.0053 / 4.5e-5 for Euler 1 /
       Euler 2 / midpoint / Gauss 4.
     - The midpoint convergence study below uses the same field.
   * - ``test/tests/test_axis_crossing.f90``
     - Six near-axis starts at two resolutions; energy excursion
       :math:`\le` 1e-3 / 1e-2, half-trace mean-energy difference
       :math:`\le` 1e-6 / 1e-4, no axis loss.
     - Not run.
   * - ``test/tests/test_orbit_invariants.f90``
     - Analytic definitions and inverse reconstruction of :math:`\bar H`,
       :math:`J_\perp`, :math:`p_\varphi`; canonical normalisation; POTATO flux
       convention. Scaled tolerances 1e-11 / 1e-10 / 1e-9.
     - Run; assertions bound unprinted residuals.
   * - ``test/tests/test_profiles.f90``
     - Flat-profile and scalar collision coefficients agree to 1e-12;
       slowing-down residence histogram matches :math:`p^2/(p^3+p_c^3)`
       within 4 %; Spitzer drag within 3 % for
       :math:`0.2 < p_c < 0.5`; Maxwellian fixed point with 8000 markers and
       4000 steps: mean :math:`p^2` within 5 %, populated-bin deviation below
       10 %.
     - Run; distribution and moment errors are printed.
   * - ``test/tests/test_fusion_source.f90`` (``test_fusion_*``, unit)
     - Bosch–Hale reactivity at 10 and 100 keV to 1e-12; fixed-seed
       :math:`\chi^2` (20 bins, 19 degrees of freedom, limit 65) for flat and
       peaked analytic profiles; the peaked source has CDF
       :math:`(s - s^2/2 - s^3/3 + s^4/4)/(5/12)`; angular and pitch moments,
       energy, Simpson power quadrature, second-order grid convergence,
       invalid inputs.
     - Not run.
   * - ``test/tests/test_fusion_birth_equilibrium.f90``
       (``test_fusion_*``, integration)
     - Flat profiles on the test VMEC equilibrium and the analytic Boozer
       chartmap: :math:`\chi^2(19) \le 65` against bin volumes from the wout
       ``vp`` (:math:`dV/ds = 4\pi^2 v_p`) or an independent chartmap
       quadrature, and total volume against ``volume_p``. A peaked VMEC case
       with three million markers must have an expected :math:`\chi^2`
       excess above 65 when the Jacobian is omitted. CDF, power and
       fixed-seed histogram are compared at 1 and 8 threads.
     - Not run.
   * - ``test/tests/test_boozer_chartmap.f90``, ``test_chartmap_scaling.f90``
     - Analytic fixture :math:`A_\vartheta = \psi_\text{tor}\rho^2`,
       :math:`h_r = 0`; :math:`\rho \leftrightarrow s` round trip
       :math:`\le` 1e-14; field, length and flux scaling.
       ``test_chartmap_rz_consistency.f90``: VMEC and chartmap boundaries
       agree within 1.5 cm (not run in the report).
     - Run; assertions bound unprinted residuals.
   * - ``test/tests/test_boozer_chartmap_roundtrip.f90``
     - VMEC → Boozer chartmap → reimport: relative :math:`B` error
       :math:`\le` 8e-5 at 50 points; 500-step midpoint orbits with
       :math:`|\Delta s| \le` 2e-9. Both paths share the equilibrium.
     - Run; both errors are printed.
   * - ``test/tests/test_spectre_sympl_volume.py``,
       ``field_can/test_magfie_spectre_analytic.f90``
     - Midpoint/Euler energy excursion and fitted secular growth below 5e-6;
       cross-scheme and short RK45 comparisons; drift derivatives against
       central differences to 1e-6.
     - Not run.
   * - ``test/tests/test_fo_boris.f90``,
       ``field_can/test_boozer_spline_derivatives.f90``
     - Full-orbit energy and adiabatic invariant; spline first and second
       derivatives against finite differences to 2e-6 and 5e-3.
     - Not run.
   * - ``tools/step_convergence.py``, ``test/python/test_step_convergence.py``
     - Paired outcomes at two step sizes: bias :math:`\le` one binomial
       standard error, flip rate :math:`\le` 3 standard errors; checked on
       synthetic ``times_lost.dat`` files.
     - Not run.
   * - ``test/python/test_simple_api.py``, ``test_orbit_macro.py``
     - ``pysimple`` and ``simple.x`` agree (shared code, not an oracle).
     - Not run.
   * - ``test/golden_record/``
     - Whole-run comparison with a reference revision (Boozer, canonical,
       Meiss, Albert, collisions, tokamak, classifier) built with the pinned
       libneo; ``compare_files.py`` rtol 1e-7 (CTest 1.5e-7), atol 1e-12.
       Regression, not an analytic oracle.
     - Not run: needs a second build.

Generate the report
-------------------

From the repository root, with NumPy, Matplotlib and netCDF4 (pytest and
SciPy for the helper tests):

.. code-block:: console

   make build-deterministic CONFIG=Fast FLAGS='-DENABLE_PYTHON_INTERFACE=OFF'
   make test CONFIG=Fast TEST='^test_verification_report$$'
   python tools/generate_verification_report.py

The report ``build/output/verification-report/verification-report.md`` sits
next to its figures, ``measurements.json``, CSV files and per-run inputs,
trajectories and logs. The budget is 300 s (``--timeout``); the output
directory must be empty (``--output-dir`` for further runs). Any missing
test, failed run, non-finite value, unresolved marker or timeout gives a
non-zero exit code. CI builds the report on dispatch with
``verification_report=true``.

Reading the figures
-------------------

*Convergence.* Midpoint endpoint error against a 256-step Gauss 2 reference
(``integmode = 5``) on the analytic TEST tokamak at normalised time 0.1; a
128-step Gauss 2 run must agree with the reference to 1 % of the smallest
midpoint error. The error is the larger absolute difference of the two
unwrapped endpoint angles, because Gauss stores the last stage radius
(:doc:`/model/integrators`). The fitted order over 8, 16, 32 and 64 steps
must lie in [1.8, 2.2]; relative energy excursion stays below 1e-3, and
:math:`J_\perp` and :math:`p_\varphi` excursions below 1e-10 (TEST units
:math:`B_0 = R_0 = \rho_0 = 1`, :math:`a = 0.5`). RK45 is excluded because
the TEST field backend in ``src/magfie.f90`` omits its curvature drift.

*Sampling uncertainty.* Three seeds of random start angles and pitch on
:math:`\rho = 0.85` of the analytic chartmap fixture, traced for 1 µs, with
nested prefixes of 32 to 256 markers. The report gives loss counts, the
binomial standard error (:doc:`/user-guide/postprocessing`), 95 % Wilson
intervals checked against SciPy, and the spread across seeds at
:math:`N = 256`. These assume independent, equally weighted markers and
exclude numerical bias and model error.
