Orbit classification
====================

SIMPLE can label collisionless orbits as regular or chaotic from a short part
of their trajectory. Use it to stop tracing regular trapped orbits early
(about ten times faster loss computations, Albert, Kasilov and Kernbichler
2020) or as a confinement metric for optimisation (Albert et al. 2023).

Method: :doc:`/model/classification`.

Classifiers
-----------

Codes per marker: 0 not classified, 1 regular, 2 chaotic (also 1 for deeply
passing markers skipped by ``contr_pp``).

:math:`J_\parallel`
   Conservation of the parallel adiabatic invariant between successive
   banana tips of trapped orbits; needs ``nturns`` tip recurrences.
Topological ("ideal orbit")
   Whether the banana-tip footprints in the Poincaré section lie on a closed
   curve; needs ``nturns`` tip recurrences.
Fractal
   Box-counting dimension of the orbit's Poincaré cut accumulated until
   ``tcut``; needs much longer tracing.

Settings
--------

===============  ============  ===================================================
Key              Default       Meaning
===============  ============  ===================================================
``tcut``         ``-1d0``      cut time in s; ``> 0`` enables the classifiers
``nturns``       8             tip recurrences for :math:`J_\parallel` and topology
``fast_class``   ``.false.``   stop trapped orbits that :math:`J_\parallel` finds
                               regular; they count as confined (not with
                               ``class_plot``)
``class_plot``   ``.false.``   also enables the classifiers; moves starts to the
                               cut and writes plot data to ``fort.*`` files
``cut_in_per``   ``0d0``       toroidal cut position as a fraction of a period
===============  ============  ===================================================

Classification cannot be combined with collisions or GPU tracing. With
classification on, SIMPLE writes ``class_parts.dat`` (index, start radius, :math:`J_\perp`, codes for
:math:`J_\parallel`, topological, fractal) and fills ``iclass`` in
``results.nc``.

Python
------

.. code-block:: python

   import pysimple

   pysimple.init("wout.nc", isw_field_type=2, trace_time=0.1)
   start = pysimple.sample_surface(128, s=0.6)
   c = pysimple.classify_fast(start, nturns=8)      # J_par and topology
   f = pysimple.classify_fractal(start, tcut=0.1)   # fractal dimension

Each result adds ``passing``, ``lost``, ``jpar``, ``topology`` and
``fractal`` to the trace results. Complete scripts: ``examples/classify_fast.py``
and ``examples/classify_fractal.py`` (:doc:`/examples`);
``examples/classification/`` reproduces the figures of Albert et al. (2023).
:doc:`/tutorials/optimisation-loop` uses ``classify_fast`` as an objective.

References
----------

* C. G. Albert, S. V. Kasilov, W. Kernbichler, J. Plasma Phys. 86,
  815860201 (2020), https://doi.org/10.1017/S0022377820000203
* C. G. Albert, R. Buchholz, S. V. Kasilov, W. Kernbichler, K. Rath,
  J. Plasma Phys. 89, 955890301 (2023), https://doi.org/10.1017/S0022377823000351
