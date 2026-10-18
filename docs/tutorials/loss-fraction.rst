Tutorial 1: Loss fraction from a VMEC equilibrium
=================================================

This tutorial computes the collisionless alpha-particle loss fraction of two
optimised reactor-scale stellarators and compares them: the
quasi-axisymmetric (QA) and the quasi-helically symmetric (QH) configuration
of Landreman and Paul (2022). Both VMEC equilibria are public simsopt test
files with a mean field strength of about 5.9 T.

You will learn how to

* prepare a run directory with a VMEC ``wout`` file and a ``simple.in``,
* choose the number of particles, tracing time and starting surface,
* read the loss fraction and the loss times,
* estimate the statistical error.

Runtime: about one minute per configuration on 32 threads.

Step 1: Input data
------------------

.. code-block:: bash

   mkdir -p tutorial1/qa tutorial1/qh
   cd tutorial1
   wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc \
       -O qa/wout.nc
   wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc \
       -O qh/wout.nc

Step 2: Input file
------------------

Write ``tutorial1/simple.in``:

.. code-block:: fortran

   &config
   netcdffile = 'wout.nc'        ! VMEC equilibrium
   isw_field_type = 2            ! Boozer canonical coordinates
   integmode = 1                 ! symplectic Euler
   trace_time = 1d-1             ! s
   sbeg = 0.6d0                  ! start surface, normalized toroidal flux
   ntestpart = 256               ! number of alpha particles
   deterministic = .true.        ! reproducible sampling
   output_results_netcdf = .true.
   /

The choices are:

``trace_time = 0.1``
   The tracing time in seconds. The slowing-down time of 3.5 MeV alphas in a
   reactor is of the order of 0.1 to 1 s. Most losses of an optimised
   configuration occur within the first 10 to 100 ms. Without collisions
   (the default), the result is the collisionless loss fraction up to this
   time.
``sbeg = 0.6``
   Particles start on the surface :math:`s = 0.6` at random positions along
   one field line (``startmode = 1``) with an isotropic pitch distribution.
   Half radius in :math:`\sqrt{s}` is :math:`s = 0.25`; a larger ``sbeg``
   probes the outer region, where losses are larger. For a volume
   distribution use ``startmode = 5`` (:doc:`/user-guide/running`).
``ntestpart = 256``
   The statistical error of a loss fraction :math:`f` from :math:`N`
   particles is :math:`\sqrt{f(1-f)/N}`. Production runs use several thousand
   particles.
``isw_field_type = 2``, ``integmode = 1``
   Canonical Boozer coordinates with the symplectic Euler integrator, the
   setting used in this tutorial (:doc:`/model/integrators`).

Copy it to both directories:

.. code-block:: bash

   cp simple.in qa/
   cp simple.in qh/

Step 3: Run
-----------

.. code-block:: bash

   (cd qa && ../../build/simple.x simple.in > simple.log)
   (cd qh && ../../build/simple.x simple.in > simple.log)

The path to ``simple.x`` assumes that ``tutorial1`` is in the repository
root. Use ``OMP_NUM_THREADS`` to set the number of threads.

Step 4: Loss fraction
---------------------

The last line of ``confined_fraction.dat`` holds the confined fractions of
passing and trapped particles at the end of the run:

.. code-block:: bash

   tail -1 qa/confined_fraction.dat
   tail -1 qh/confined_fraction.dat

::

   0.10000000000000001  0.70703125000000000  0.17187500000000000  256
   0.10000000000000001  0.67968750000000000  0.32031250000000000  256

For QA the confined fraction is :math:`0.707 + 0.172 = 0.879`, so 12.1 % of
the particles are lost within 0.1 s. In QH no particle is lost.

``orbit_exit_code.dat`` gives the reason why each orbit ended. Count the codes
in the second column:

.. code-block:: bash

   awk '{print $2}' qa/orbit_exit_code.dat | sort | uniq -c

::

        95 0
        31 1
       130 3

Code 0 means the particle was traced to ``trace_time``, code 1 that it crossed
the last closed flux surface, and code 3 that it was not traced. Code 3
marks strongly passing particles: by default, particles with trapping
parameter below ``contr_pp = -1`` are counted as confined without tracing,
because they stay close to their flux surface. Set ``contr_pp = -1d10`` to
trace all particles. Codes of 100 and above mark numerically unresolved
orbits, which are counted separately and reported in
``unresolved_fraction.dat`` (:doc:`/user-guide/outputs`).

Step 5: Plot the confined fraction over time
--------------------------------------------

Save as ``tutorial1/plot_confined.py`` and run it in ``tutorial1`` with
``python plot_confined.py`` (``matplotlib`` is installed by
``setup-venv.sh``):

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np

   for case in ("qa", "qh"):
       t, passing, trapped, resolved = np.loadtxt(
           f"{case}/confined_fraction.dat", unpack=True)
       loss = 1.0 - (passing + trapped)
       n = resolved[-1]
       error = np.sqrt(loss[-1] * (1.0 - loss[-1]) / n)
       print(f"{case}: loss fraction {loss[-1]:.3f} +- {error:.3f}")
       plt.semilogx(t, loss, label=case.upper())

   plt.xlabel("time (s)")
   plt.ylabel("loss fraction")
   plt.legend()
   plt.savefig("loss_fraction.png", dpi=150)

::

   qa: loss fraction 0.121 +- 0.020
   qh: loss fraction 0.000 +- 0.000

The plot shows when the losses happen. Zero losses from 256 independent
markers give a one-sided 95 % binomial upper bound of
:math:`1 - 0.05^{1/N} \simeq 3/N = 1.2\,\%`, rather than evidence that the
loss probability vanishes. The plug-in standard error printed above is zero
in this case and does not express that uncertainty.

Step 6: Loss times of individual particles
------------------------------------------

``times_lost.dat`` has one line per particle: index, loss time in seconds,
trapping parameter, starting :math:`s`, perpendicular invariant, and the five
phase-space coordinates at the end. Confined particles have the loss time
``trace_time``; skipped particles have ``-1``. For example, the five earliest
losses in QA:

.. code-block:: bash

   awk '$2 > 0 && $2 < 0.1' qa/times_lost.dat | LC_ALL=C sort -g -k2 | head -5

The trapping parameter (third column) is negative for passing and between 0
and 1 for trapped particles. In this QA run 28 of the 31 lost particles are
trapped. The first losses occur after about 20 microseconds.

What next
---------

* Change ``sbeg``, ``ntestpart`` and ``trace_time`` and observe the loss
  fraction and its error.
* Add collisions with ``swcoll = .true.`` and background profiles
  (:doc:`/user-guide/running`).
* :doc:`wall-heat-loads` uses the QA result file ``results.nc`` to compute a
  heat-load pattern.
