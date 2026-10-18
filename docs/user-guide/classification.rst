Orbit classification
====================

Without collisions, a guiding-centre orbit on a regular (non-chaotic) orbit
stays confined for all times. SIMPLE can classify orbits as regular or chaotic
from a short part of their trajectory. This has two uses:

* **Acceleration** of loss computations: regular trapped orbits are counted as
  confined and their tracing stops early (``fast_class``). Albert, Kasilov and
  Kernbichler (2020, JPP) report a speed-up of about one order of magnitude.
* **Confinement metrics** for optimisation: the fraction of chaotic orbits,
  and the extent of regular regions that act as barriers, measure the alpha
  confinement with much shorter tracing (Albert et al. 2023).

The methods are described in :doc:`/model/classification`.

Classifiers
-----------

Three classifiers are available. Each gives 0 (not classified), 1 (regular)
or 2 (chaotic) per particle.

:math:`J_\parallel` classifier
   Uses the conservation of the parallel adiabatic invariant at the banana
   tips of trapped orbits. Needs ``nturns`` tip recurrences.
Topological classifier
   Tests whether the footprints of the banana tips in the Poincaré section
   lie on a closed curve, with recurrence and monotonicity conditions
   ("ideal" orbits). Needs ``nturns`` tip recurrences.
Fractal classifier
   Estimates the Minkowski (box-counting) dimension of the Poincaré cut of the
   orbit, accumulated until the cut time ``tcut``. Needs much longer tracing.

Settings
--------

.. code-block:: fortran

   &config
   nturns = 8            ! tip recurrences for the J_par and topological test
   fast_class = .true.   ! stop regular trapped orbits after classification
   tcut = 1d-1           ! fractal classification at this time (s); -1 = off
   class_plot = .false.  ! write Poincare plot data (fort.* files)
   cut_in_per = 0d0      ! position of the cut in the field period, [0, 1]
   /

``simple.x`` runs the classifiers when ``tcut > 0`` or
``class_plot = .true.``; ``fast_class`` then stops the tracing of regular
trapped orbits. The classifiers use the symplectic integrators
(``integmode > 0``) and cannot be combined with collisions. With
``tcut > 0`` or ``class_plot = .true.``, SIMPLE writes
``class_parts.dat``: index, start radius, perpendicular invariant and the
three classes (:math:`J_\parallel`, topological, fractal).

Python
------

.. code-block:: python

   import pysimple

   pysimple.init("wout.nc", isw_field_type=2, trace_time=0.1)
   start = pysimple.sample_surface(128, s=0.6)
   c = pysimple.classify_fast(start, nturns=8)      # J_par and topology
   f = pysimple.classify_fractal(start, tcut=0.1)   # fractal dimension

The result dictionaries contain ``passing``, ``lost``, ``jpar``,
``topology`` and ``fractal`` in addition to the trace results.
:doc:`/tutorials/optimisation-loop` uses ``classify_fast`` as a confinement
metric. ``examples/classify_fast.py`` and ``examples/classify_fractal.py``
(:doc:`/examples`) are complete scripts, and
``examples/classification/`` reproduces the figures of Albert et al. (2023).
