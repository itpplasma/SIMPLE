Orbit classification
====================

This page defines the footprints and the three classifiers that label an
orbit regular or chaotic from a short trace, and the class codes written to
the outputs. Usage: :doc:`/user-guide/classification`.

Without collisions only chaotic orbits are lost late, so a regular orbit
can be stopped early (Albert, Kasilov and Kernbichler 2020, fractal
classifier; Albert et al. 2023, :math:`J_\parallel` and topological
classifiers and metrics). Code: ``src/classification.f90``,
``src/check_orbit_type.f90``, ``src/cut_detector.f90``. Classification
excludes collisions and full-orbit tracing.

Footprints
----------

The classifier keeps a stencil of the last ``nplagr`` :math:`= 4`
microstep states :math:`(r,\vartheta,\varphi,v/v_0,\lambda)` together with the
running integral

.. math::

   J_\parallel = \int \lambda^2\, d\tau ,

which, at constant speed, is proportional to the parallel adiabatic invariant
:math:`\oint m v_\parallel\, dl` over one bounce.

* **Banana tips**: when :math:`\lambda` changes sign from negative to
  positive, the state is interpolated to :math:`\lambda = 0` by a Lagrange
  polynomial in :math:`\lambda` through the stencil. The footprint is
  :math:`(r,\vartheta)` at the tip, together with :math:`J_\parallel` since
  the previous tip; :math:`J_\parallel` is then reset.
* **Field-period planes**: when :math:`\varphi` crosses a multiple of
  :math:`2\pi/N_\text{fp}`, the state is interpolated to the plane in the same
  way; the footprint is :math:`(r,\vartheta)`.

Tip footprints feed all three classifiers; plane footprints only the fractal
one.

Class codes
-----------

``class_parts.dat`` and ``results.nc`` store three codes per marker
(:doc:`/user-guide/outputs`): :math:`J_\parallel`, topological, fractal.
Each is 0 (unclassified), 1 (regular) or 2 (chaotic). Passing markers skipped
by ``notrace_passing`` or ``contr_pp`` get 1 in all three.

:math:`J_\parallel` and topological classifiers
-----------------------------------------------

``check_orbit_type`` processes the tip footprints in order. A *recurrence*
occurs when the angle of a new footprint lies between the angles of the
first two footprints (angles taken modulo :math:`2\pi` around the first).
If the first recurrence comes after at most 10 footprints, both codes are 2.
Otherwise, once ``nturns`` (default 8) recurrences are collected:

* **Topological**: the numbers of footprints between consecutive recurrences
  must be consistent with one closed curve. For every multiplicity
  :math:`m \le \texttt{nturns}/2`, each interval of :math:`m` recurrences must
  contain :math:`N_m` or :math:`N_m - 1` footprints, where :math:`N_m` is the
  first such interval (``check_recurs``). In addition, the footprints of each
  recurrence, sorted by angle, must keep their order in the following
  recurrences. Both conditions hold: code 1 (ideal); otherwise 2.
* :math:`J_\parallel`: code 1 if :math:`|J_\parallel^{(k)} - J_\parallel^{(2)}|
  \le 15` for every tip :math:`k` up to the last recurrence, otherwise 2.
  :math:`J_\parallel` is in units of :math:`\tau` (cm), so the tolerance is
  absolute.

With ``fast_class = .true.`` and without ``class_plot``, a trapped marker is
stopped and counted confined once its :math:`J_\parallel` code is 1.

Fractal classifier
------------------

At the microstep ``ntcut`` (from ``tcut``) the :math:`N` collected footprints
of each kind are box-counted (``fract_dimension``) on their bounding
rectangle in :math:`(r,\vartheta)` with a :math:`2^n \times 2^n` grid,
:math:`n = \lfloor\log_4 N\rfloor`, so that the number of boxes is about
:math:`N`. The occupied fraction

.. math::

   \phi = N_\text{boxes}/4^n

is :math:`\sim 2^{-n}` for footprints on a curve and of order one for an
area; :math:`\phi > 0.2` (box-counting dimension
:math:`d > 2 + \log 0.2/\log 2^n`) is chaotic. Tip footprints set the fractal
code; chaotic plane footprints only mark the orbit as not regular. This
classifier needs far longer traces but also covers passing orbits.

Albert et al. (2023) derive two metrics from the codes: the fraction of
chaotic trapped orbits, and the radial extent of regular regions in
:math:`J_\perp`, which bound chaotic orbits starting further inside.

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
