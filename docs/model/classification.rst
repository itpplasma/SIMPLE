Orbit classification
====================

Without collisions, a guiding-centre orbit that lies on an invariant torus
(a regular orbit) is confined for all times if the torus lies inside the
plasma. Only chaotic orbits can be lost after a long time, and their loss time
is unpredictable. SIMPLE classifies orbits as regular or chaotic from a short
piece of trajectory. Albert, Kasilov and Kernbichler (2020, JPP) used a
fractal-dimension classifier to stop regular orbits early; Albert et al.
(2023) introduced two faster classifiers based on the parallel adiabatic
invariant and on the topology of the orbit footprints, and confinement
metrics derived from them. The implementation is in ``src/classification.f90``,
``src/check_orbit_type.f90`` and ``src/cut_detector.f90``.

Footprints in Poincaré sections
-------------------------------

Two kinds of Poincaré sections are used.

* **Banana tips** of trapped orbits: the points where :math:`v_\parallel`
  changes sign from negative to positive. Each tip is located by Lagrange
  interpolation over a stencil of integration points around the sign change.
  Its footprint is the radial position and an angle at the tip.
* **Field-period planes** :math:`\varphi = \text{const}` (offset
  ``cut_in_per`` within the period) for passing orbits.

Between two consecutive tips the code accumulates the parallel adiabatic
invariant

.. math::

   J_\parallel \propto \oint v_\parallel^2 \, d\tau ,

which is proportional to :math:`\oint m v_\parallel \, dl` along the orbit.

:math:`J_\parallel` classifier
------------------------------

For a regular trapped orbit, the tip footprints lie on a closed curve, and
:math:`J_\parallel` is a smooth function along this curve, because it is an
adiabatic invariant of the bounce motion. The classifier sorts the footprints
of ``nturns`` recurrences along the curve and tests whether
:math:`J_\parallel` varies smoothly within a tolerance. If so, the orbit is
regular (class 1), otherwise stochastic (class 2). Orbits that do not complete
enough recurrences remain unclassified (class 0).

Topological classifier
----------------------

The same footprints are tested for the topology of an "ideal" regular orbit:
the footprints must recur on a single closed line, and their order along the
line must advance monotonically from one recurrence to the next. Orbits that
pass both tests are ideal (class 1), others non-ideal (class 2).

Both classifiers need only a few bounce periods (``nturns`` recurrences,
default 8). With ``fast_class = .true.``, the tracing of a trapped orbit stops
when the :math:`J_\parallel` classifier finds it regular, and the orbit is
counted as confined.

Fractal classifier
------------------

At the cut time ``tcut`` the collected footprints (tips for trapped orbits,
field-period crossings for passing orbits) are analysed with a box-counting
method that estimates their Minkowski dimension: the footprints of a regular
orbit fill a curve (dimension 1), those of a chaotic orbit fill an area
(dimension 2). The measure is compared with a fixed threshold (0.2 in
``src/classification.f90``). This classifier needs much longer tracing than
the other two, but it also applies to passing orbits.

Confinement metrics
-------------------

The classes give metrics for collisionless losses that need far less tracing
time than a direct loss computation (Albert et al. 2023):

* the fraction of chaotic trapped orbits on a surface or in the volume,
* the radial extent of regular regions in the space of the perpendicular
  invariant, which act as barriers for chaotic orbits that start further
  inside.

Because pitch-angle scattering becomes important only after the alphas have
lost a large part of their energy, these metrics remain informative with
collisions.

References
----------

* C. G. Albert, S. V. Kasilov, W. Kernbichler, Accelerated methods for direct
  computation of fusion alpha particle losses within stellarator
  optimization, J. Plasma Phys. 86, 815860201 (2020),
  https://doi.org/10.1017/S0022377820000203
* C. G. Albert, R. Buchholz, S. V. Kasilov, W. Kernbichler, K. Rath, Alpha
  particle confinement metrics based on orbit classification in
  stellarators, J. Plasma Phys. 89, 955890301 (2023),
  https://doi.org/10.1017/S0022377823000351
