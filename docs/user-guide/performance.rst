Performance
===========

This page lists what sets the run time of SIMPLE and how to reduce it. A run
has two phases, both timed in the log: field setup (splines, canonical or
Boozer transformation, start surface) and tracing. Setup dominates for a few
hundred markers, tracing for thousands. Measured scaling and profile:
:doc:`/performance/efss-hpc/index`.

Build and threads
-----------------

* Production: ``make CONFIG=Fast`` (native CPU tuning). Use
  ``-DSIMPLE_DETERMINISTIC_FP=ON`` (no fast-math) only for regression
  comparisons.
* Markers are traced in parallel with OpenMP. Set ``OMP_NUM_THREADS`` to the
  physical cores, ``OMP_PLACES=cores`` and ``OMP_PROC_BIND=spread``. Tracing
  scales nearly linearly up to 16 threads; setup is mostly serial.
* One run uses one node. Spread large ensembles over independent jobs
  (equilibria, seeds, or ``batch_size`` batches); cluster scripts are in
  :doc:`/hpc/index`.
* For many markers in Boozer coordinates without collisions, walls or
  classification, use the CUDA backend (:doc:`gpu`).

Physics settings
----------------

* Cost is proportional to ``ntestpart`` × ``trace_time`` for confined
  markers; lost markers stop early.
* ``contr_pp`` (default :math:`-1`) counts strongly passing markers as
  confined without tracing. In :doc:`/tutorials/loss-fraction` this skips 130
  of 256 markers; ``contr_pp = -1d10`` traces them.
* ``fast_class = .true.`` with ``tcut > 0`` stops regular trapped orbits after
  classification (:doc:`classification`).
* ``macrostep_time_grid = 'log'`` resolves early losses in time without
  adding microsteps.

Numerical settings
------------------

* ``npoiper2`` (default 256) sets the step size. Check the loss fraction for
  convergence by halving the step; symplecticity does not guarantee stability
  at large steps.
* Higher-order integrators (``integmode`` 3 to 7, 15) cost more per step;
  symplectic Euler (``integmode = 1``) is usually fastest for loss fractions.
* ``multharm``, ``ns_s`` and ``ns_tp`` set the VMEC spline resolution; lower
  values (for example ``multharm = 3``) speed up setup for scans.
* ``canonical_grid_nr``, ``canonical_grid_ntheta``, ``canonical_grid_nphi``
  set the Meiss and Albert transformation grid; change one at a time in
  resolution studies.

Guidance per key: :doc:`/reference/numerical-settings`.

.. toctree::
   :maxdepth: 1

   /performance/efss-hpc/index
