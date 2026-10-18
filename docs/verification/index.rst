Verification and benchmarks
===========================

.. placeholder-page

.. note::

   This page will hold a verification report generated from the test suite
   and the benchmark cases. Until it is added, the page summarises the
   published benchmarks and the tests that run with every change.

Published benchmarks
--------------------

The methods of SIMPLE have been verified and validated in peer-reviewed
publications:

* **Symplectic integrators.** Albert, Kasilov and Kernbichler (2020, J.
  Comput. Phys.) study guiding-centre orbits in an axisymmetric tokamak and a
  realistic three-dimensional stellarator field. They show the better
  long-term properties of the symplectic methods compared with an adaptive
  Runge-Kutta scheme. For alpha-particle loss statistics over the
  slowing-down time, symplectic Euler is more than three times faster at the
  same statistical accuracy, with linear scaling in the number of threads.
* **Loss fractions.** Albert, Kasilov and Kernbichler (2020, J. Plasma Phys.)
  compute collisionless alpha losses for an ensemble of 1000 orbits in a
  quasi-isodynamic, a quasi-helical and a quasi-axisymmetric configuration and
  combine the symplectic integrators with early classification. They obtain
  reliable results at about one order of magnitude lower cost than direct
  integration with adaptive Runge-Kutta methods.
* **Classification metrics.** Albert et al. (2023, J. Plasma Phys.) show that
  the :math:`J_\parallel` and topological classifiers need at least one order
  of magnitude less tracing time than the fractal classifier, and compare the
  resulting metrics with collisional loss computations for two
  quasi-isodynamic, one quasi-helical and two quasi-axisymmetric
  configurations.
* **Wall heat loads.** Albert et al. (2026, Nucl. Mater. Energy) compute
  alpha-particle wall heat loads in stellarator reactors with the STL wall
  intersection and compare them with the BEAMS3D code.

The full references are in :doc:`/reference/publications`.

Automated tests
---------------

Every pull request runs the test suite on GitHub Actions. ``make test``
runs it locally (:doc:`/development/index`). It contains

* unit tests of splines, coordinate transformations, field evaluation,
  canonical constructions and the STL intersection,
* integration tests of the integrators, including symplecticity and
  invariant-conservation checks on analytic fields,
* Python tests that compare ``pysimple`` with ``simple.x`` and test the
  examples,
* golden-record tests, which compare the output of the current version with
  a reference version for a set of cases (``test/golden_record``,
  ``make test-golden-main``).

The review policy (:doc:`/review-policy`) defines which evidence
a change must provide, including independent comparisons for changes of
physics and output semantics.
