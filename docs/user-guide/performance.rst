Performance tips
================

A SIMPLE run has two phases: the field setup (VMEC splines, Boozer or
canonical transformation, starting surface) and the particle tracing. The log
prints the time of each phase. For a few hundred particles the setup can
dominate; for thousands of particles the tracing does.

Build
-----

* Use ``make CONFIG=Fast`` (native CPU optimisation) for production runs.
* ``SIMPLE_DETERMINISTIC_FP=ON`` disables fast-math. Use it for regression
  comparisons, not for production.

Parallel runs
-------------

* Set ``OMP_NUM_THREADS`` to the number of physical cores, and
  ``OMP_PROC_BIND=spread``, ``OMP_PLACES=cores`` on multi-socket nodes.
* For large ensembles, run independent cases (equilibria, seeds or batches
  with ``batch_size``) as separate jobs. Example Slurm scripts are in
  ``RUN/``.
* For many particles in Boozer coordinates without collisions, walls or
  classification, use the CUDA backend (:doc:`gpu`).

Physics settings that set the cost
----------------------------------

* The cost is proportional to ``ntestpart`` times ``trace_time`` for confined
  particles. Lost particles stop early.
* ``contr_pp`` (default :math:`-1`) skips strongly passing particles, which
  are well confined. In :doc:`/tutorials/loss-fraction` 130 of 256 particles
  are skipped this way; tracing them (``contr_pp = -1d10``) roughly doubles
  the number of orbits traced to the end.
* ``fast_class = .true.`` with ``tcut > 0`` stops regular trapped orbits after
  classification (:doc:`classification`).
* A logarithmic macrostep grid (``macrostep_time_grid = 'log'``) resolves early
  losses without increasing the number of microsteps.

Numerical settings
------------------

* ``npoiper2`` sets the step size. Symplectic Euler can give bounded long-time
  energy errors at suitable steps, but symplecticity alone does not guarantee
  stability at large steps. Check
  the loss fraction for convergence by halving the step. ``DOC/config.md``
  gives guidance (:doc:`/reference/design-notes`).
* Higher-order integrators (``integmode`` 3 to 7, 15) are more accurate per
  step but cost more per step. Symplectic Euler is usually the fastest
  choice for loss fractions.
* ``multharm``, ``ns_s`` and ``ns_tp`` set the spline resolution of the VMEC
  field. Lower values speed up the setup for quick scans (for example
  ``multharm = 3``).
* ``canonical_grid_nr``, ``canonical_grid_ntheta`` and
  ``canonical_grid_nphi`` set the grid of the Meiss and Albert
  transformation. Change one at a time in resolution studies.

The :doc:`/performance/efss-hpc/index` gives measured OpenMP scaling,
profiles and commands for reproducing the comparison.
