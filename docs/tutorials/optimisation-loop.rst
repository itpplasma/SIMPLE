Tutorial 3: SIMPLE inside an optimisation loop
==============================================

Stellarator optimisation evaluates many candidate equilibria. This tutorial
wraps a SIMPLE loss computation into a Python function that takes the path of
a VMEC ``wout`` file and returns the loss fraction with its statistical error.
Such a function can be called from any optimiser, for example after each
equilibrium solve in simsopt.

You will learn how to

* configure and run SIMPLE from Python with ``pysimple``,
* write a reproducible loss objective with common random numbers,
* judge whether a change of the objective is significant,
* use the fast orbit classification as a cheaper proxy.

Prerequisites: ``pysimple`` installed with ``setup-venv.sh``
(:doc:`/getting-started/installation`) and the two equilibria of
:doc:`loss-fraction` in ``tutorial1/qa`` and ``tutorial1/qh``. Runtime: about
one minute.

The objective function
----------------------

.. literalinclude:: scripts/loss_objective.py
   :language: python
   :caption: docs/tutorials/scripts/loss_objective.py

``pysimple.init`` resets all ``&config`` values to their defaults, then
applies the keyword arguments. Its first argument is the equilibrium file. It
is used as ``netcdffile``, ``field_input`` and ``coord_input``. Every call
loads a new equilibrium, so the same Python process can evaluate any number of
candidates.

``pysimple.sample_surface`` places the start points on one field line of the
surface :math:`s` = ``S_START``. With ``deterministic=True`` the random
generator is reseeded in every ``init``. Two calls with the same equilibrium
therefore use the same start points, and two similar equilibria use nearly the
same ones. This technique of common random numbers removes most of the sampling
noise from the *difference* of two objective values, which is what an
optimiser compares.

``trace_parallel`` traces all particles with OpenMP and returns the loss time
of each one: ``trace_time`` if it stayed confined, ``NaN`` if the orbit could
not be resolved numerically, ``-1`` if it was skipped as a deeply passing
particle. The function counts unresolved orbits neither as lost nor as
confined.

Run
---

.. code-block:: bash

   source .venv/bin/activate
   mkdir -p tutorial3
   cd tutorial3
   python ../docs/tutorials/scripts/loss_objective.py \
       ../tutorial1/qa/wout.nc ../tutorial1/qh/wout.nc

The script prints the SIMPLE log and then::

   QA, B x 0.50: loss fraction 0.344 +- 0.042  (11 s)
   QA, B x 0.75: loss fraction 0.242 +- 0.038  (12 s)
   QA, B x 1.00: loss fraction 0.141 +- 0.031  (12 s)
   QH, B x 1.00: loss fraction 0.000 +- 0.000

The first three lines scan a single design parameter, the field strength,
through ``vmec_B_scale``, which multiplies the VMEC field. A weaker field
means larger orbit widths and more losses. The last line evaluates the QH
equilibrium with the same settings. Each evaluation takes about ten seconds,
of which a few seconds are the field setup.

Is a change significant?
------------------------

The standard error :math:`\sqrt{f(1-f)/N}` holds for one evaluation. For a
difference between two equilibria with common random numbers the error is
smaller, because the same particles are compared. In an optimisation:

* use a fixed seed (``deterministic=True``) and the same ``N`` for all
  evaluations,
* increase ``N`` as the optimisation converges and the differences shrink,
* do not trust differences much below the standard error of one evaluation
  unless you have checked their noise by repeating with another seed
  (``ran_seed``).

A shorter ``trace_time`` reduces the cost in proportion. Losses of good
configurations happen mainly in the first milliseconds (:doc:`loss-fraction`),
so 10 ms are a common choice inside an optimisation, followed by longer runs
for the final candidates.

A cheaper proxy: orbit classification
-------------------------------------

The classifiers of Albert et al. (2023) decide after a few bounce periods
whether a trapped orbit is regular (confined without collisions) or chaotic
(:doc:`/model/classification`). ``classify_fast`` stops after ``nturns``
recurrences of the banana tip:

.. literalinclude:: scripts/classify_proxy.py
   :language: python
   :caption: docs/tutorials/scripts/classify_proxy.py

.. code-block:: bash

   python ../docs/tutorials/scripts/classify_proxy.py ../tutorial1/qa/wout.nc

::

   classified in 4 s
   trapped: 43, J_par regular: 2, chaotic: 17, unclassified: 24
   lost during classification: 18

The fraction of chaotic trapped orbits is a metric for losses that needs a
fraction of the tracing time. Orbits that are lost before they complete
``nturns`` recurrences remain unclassified.

Embedding in an optimiser
-------------------------

The function ``loss_fraction`` only needs the path of a ``wout`` file. In an
optimisation framework, call it after each equilibrium solve with the
``wout`` file written by the solver and add the result, scaled by a weight, to
the objective. Keep three points in mind:

* The objective is noisy and not differentiable with respect to the
  equilibrium. Use derivative-free optimisers or smooth proxies (for example
  quasi-symmetry residuals) for the main optimisation, and the loss fraction
  as a constraint or as a term with a modest weight.
* SIMPLE uses all OpenMP threads of the process. Set ``OMP_NUM_THREADS`` when
  the optimiser runs several evaluations in parallel.
* An error inside the Fortran library can stop the whole Python process
  (Fortran ``stop``). For long unattended optimisations, run each evaluation
  in a separate process, for example with ``multiprocessing`` or a
  ``subprocess`` call of a small script, so that one failed equilibrium does
  not end the optimisation.

For large particle numbers inside optimisation loops, the CUDA backend
computes the segmented loss objective of Landreman et al. on a GPU
(:doc:`/user-guide/gpu`).
