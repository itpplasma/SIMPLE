Verification
============

This section lists the published benchmarks of SIMPLE's methods and what the
automated tests check. :doc:`/verification` maps each numerical claim to its
test, oracle and tolerance and explains the generated report.

Published benchmarks
--------------------

* **Symplectic integrators.** Albert, Kasilov and Kernbichler (2020, J.
  Comput. Phys.) compare the symplectic schemes with adaptive Runge–Kutta for
  guiding-centre orbits in an axisymmetric tokamak and a three-dimensional
  stellarator field. For alpha loss statistics over the slowing-down time,
  symplectic Euler is more than three times faster at equal statistical
  accuracy.
* **Loss fractions.** Albert, Kasilov and Kernbichler (2020, J. Plasma
  Phys.) compute collisionless alpha losses of 1000 markers in
  quasi-isodynamic, quasi-helical and quasi-axisymmetric configurations; with
  early classification the cost is about one order of magnitude below
  adaptive Runge–Kutta tracing.
* **Classification metrics.** Albert et al. (2023, J. Plasma Phys.) show that
  the :math:`J_\parallel` and topological classifiers need at least ten times
  less tracing than the fractal classifier, and compare the metrics with
  collisional losses for two quasi-isodynamic, one quasi-helical and two
  quasi-axisymmetric configurations.
* **Wall heat loads.** Albert et al. (2026, Nucl. Mater. Energy) compute
  alpha wall heat loads in stellarator reactors with the STL wall
  intersection and compare with BEAMS3D (references:
  :doc:`/reference/publications`).

Automated tests
---------------

``make test`` runs all tests except the golden-record regression tests,
which ``make test-regression`` runs; labels and oracle rules are in
:doc:`/development/testing`. The suites cover

splines, coordinate maps, canonical constructions and the STL intersection
against analytic fixtures; per-scheme energy bounds, axis crossing and
invariant reconstruction on an analytic tokamak; collisions against
slowing-down, Spitzer-drag and Maxwellian oracles; fusion-birth sampling
against analytic sources; ``pysimple`` against ``simple.x``; and whole-run
golden records (``test/golden_record``).

.. toctree::
   :maxdepth: 1

   /verification
