Tutorial 3: SIMPLE inside an optimisation loop
==============================================

This tutorial wraps a SIMPLE loss computation into a Python function that
takes a VMEC ``wout`` file and returns the loss fraction with its error, ready
to be called by an optimiser after each equilibrium solve. It shows how to
make the objective reproducible with common random numbers, when a change is
significant, and how orbit classification serves as a cheaper proxy.

Prerequisites: ``pysimple`` from ``setup-venv.sh``
(:doc:`/getting-started/installation`) and the equilibria of
:doc:`loss-fraction` in ``tutorial1/qa`` and ``tutorial1/qh``. Runtime: about
one minute.

The objective function
----------------------

.. literalinclude:: scripts/loss_objective.py
   :language: python
   :caption: docs/tutorials/scripts/loss_objective.py

``pysimple.init`` resets all ``&config`` keys to their defaults, applies the
keyword arguments and loads the equilibrium given as first argument (used as
``netcdffile``, ``field_input`` and ``coord_input``). One Python process can
therefore evaluate any number of candidates.

``trace_parallel`` traces all markers with OpenMP and returns per marker the
loss time, ``trace_time`` if confined, ``NaN`` if numerically unresolved and
``-1`` if skipped as strongly passing. The function counts unresolved orbits
neither as lost nor as confined.

``S_START = 0.6`` differs from the default ``sbeg = 0.5`` on purpose
(limitation in :doc:`/user-guide/python`); keep the two different when
changing settings.

Common random numbers
---------------------

``sample_surface`` places the markers on one field line of the surface
``S_START``. With ``deterministic=True`` the generator is reseeded in every
``init``, so two calls on the same equilibrium use the same start points and
two similar equilibria nearly the same ones. An optimiser compares objective
values; for two candidates with loss indicators :math:`X_i, Y_i \in \{0,1\}`
on the same markers,

.. math::

   \operatorname{Var}(\hat f_X - \hat f_Y)
   = \frac{\operatorname{Var}X + \operatorname{Var}Y - 2\operatorname{Cov}(X,Y)}{N},

which is smaller than for independent samples whenever the outcomes are
positively correlated. Estimate it from the paired per-marker outcomes, not by
adding the two single-run errors.

Run
---

.. code-block:: bash

   source .venv/bin/activate
   mkdir -p tutorial3
   cd tutorial3
   python ../docs/tutorials/scripts/loss_objective.py \
       ../tutorial1/qa/wout.nc ../tutorial1/qh/wout.nc

After the SIMPLE log the script prints::

   QA, B x 0.50: loss fraction 0.344 +- 0.042  (11 s)
   QA, B x 0.75: loss fraction 0.242 +- 0.038  (12 s)
   QA, B x 1.00: loss fraction 0.141 +- 0.031  (12 s)
   QH, B x 1.00: loss fraction 0.000 +- 0.000

The first three lines scan the field strength through ``vmec_B_scale``, which
multiplies the VMEC field: a weaker field means wider orbits and more losses.
A few seconds of each evaluation are field setup. For the zero-loss QH line
the printed error is not an upper bound; use :math:`3/N`
(:doc:`loss-fraction`).

In an optimisation:

* keep the seed (``deterministic=True``) and ``N`` fixed across evaluations;
* increase ``N`` as the optimisation converges and differences shrink;
* do not trust differences below the single-evaluation error unless repeated
  with another ``ran_seed``;
* use ``trace_time = 1d-2`` (cost scales with it; most losses of good
  configurations occur in the first milliseconds) and longer runs for the
  final candidates.

A cheaper proxy: orbit classification
-------------------------------------

The classifiers of Albert et al. (2023) decide after a few bounce periods
whether a trapped orbit is regular or chaotic (:doc:`/model/classification`).
``classify_fast`` stops after ``nturns`` recurrences of the banana tip.
Regular motion implies confinement only when the invariant torus lies inside
the loss boundary.

.. literalinclude:: scripts/classify_proxy.py
   :language: python
   :caption: docs/tutorials/scripts/classify_proxy.py

.. code-block:: bash

   python ../docs/tutorials/scripts/classify_proxy.py ../tutorial1/qa/wout.nc

::

   classified in 4 s
   trapped: 43, J_par regular: 2, chaotic: 17, unclassified: 24
   lost during classification: 18

The fraction of chaotic trapped orbits is a loss metric at a fraction of the
tracing cost. Orbits lost before ``nturns`` recurrences stay unclassified.

Embedding in an optimiser
-------------------------

Call ``loss_fraction`` with the ``wout`` file written by the equilibrium
solver and add the result with a weight to the objective.

* The objective is noisy and not differentiable in the equilibrium. Use
  derivative-free optimisers, or smooth proxies (for example quasi-symmetry
  residuals) for the main optimisation with the loss fraction as a
  constraint or modestly weighted term.
* SIMPLE uses all OpenMP threads of the process; set ``OMP_NUM_THREADS`` when
  the optimiser runs evaluations in parallel.
* A Fortran ``stop`` inside the library ends the Python process. For long
  unattended runs, evaluate each candidate in a separate process
  (``multiprocessing`` or ``subprocess``).

For large marker numbers, the CUDA backend computes the segmented loss
objective of Landreman et al. on a GPU (:doc:`/user-guide/gpu`).
