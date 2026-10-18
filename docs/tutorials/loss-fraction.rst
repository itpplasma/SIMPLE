Tutorial 1: Loss fraction from a VMEC equilibrium
=================================================

This tutorial computes the collisionless alpha-particle loss fraction of the
quasi-axisymmetric (QA) and quasi-helically symmetric (QH) reactor-scale
configurations of Landreman and Paul (2022), and its statistical error. Both
VMEC equilibria are public simsopt test files with a mean field of about
5.9 T. Runtime: one to two minutes per configuration on 16 threads.

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

Write ``tutorial1/simple.in`` and copy it to both directories
(``cp simple.in qa/ && cp simple.in qh/``):

.. code-block:: fortran

   &config
   netcdffile = 'wout.nc'        ! VMEC equilibrium
   isw_field_type = 2            ! Boozer canonical coordinates
   integmode = 1                 ! symplectic Euler
   trace_time = 1d-1             ! s
   sbeg = 0.6d0                  ! start surface, normalised toroidal flux
   ntestpart = 256               ! number of alpha particles
   deterministic = .true.        ! reproducible sampling
   output_results_netcdf = .true.
   /

``trace_time = 0.1``
   Without collisions (the default) the result is the collisionless loss
   fraction up to this time. Losses of an optimised configuration occur
   mostly within the first 10 to 100 ms; the slowing-down time of 3.5 MeV
   alphas in a reactor is 0.1 to 1 s.
``sbeg = 0.6``
   Markers start on the surface :math:`s = 0.6` along one field line
   (``startmode = 1``, the default) with isotropic pitch. Larger ``sbeg``
   probes the outer region, where losses are larger. Other start
   distributions: :doc:`/user-guide/running`.
``ntestpart = 256``
   Sets the statistical error (below). Production runs use several thousand
   markers.

Step 3: Run
-----------

.. code-block:: bash

   (cd qa && ../../build/simple.x simple.in > simple.log)
   (cd qh && ../../build/simple.x simple.in > simple.log)

The relative path assumes ``tutorial1`` is in the repository root.
``OMP_NUM_THREADS`` sets the number of threads.

Step 4: Loss fraction
---------------------

The last line of ``confined_fraction.dat`` holds time, the confined fractions
of passing and trapped markers, and the number of numerically resolved
markers:

.. code-block:: bash

   tail -n 1 qa/confined_fraction.dat
   tail -n 1 qh/confined_fraction.dat

::

   0.10000000000000001  0.70703125000000000  0.17187500000000000  256
   0.10000000000000001  0.67968750000000000  0.32031250000000000  256

In QA the confined fraction is :math:`0.707 + 0.172 = 0.879`: 12.1 % of the
markers are lost within 0.1 s. In QH none is lost.

``orbit_exit_code.dat`` gives the reason each orbit ended (second column;
all codes in :doc:`/user-guide/outputs`):

.. code-block:: bash

   awk '{print $2}' qa/orbit_exit_code.dat | sort | uniq -c

::

        95 0
        31 1
       130 3

Code 0 is traced to ``trace_time``, code 1 crossed the last closed flux
surface, code 3 was not traced: markers with trapping parameter below
``contr_pp`` (default :math:`-1`) are strongly passing and counted as confined
without tracing. Set ``contr_pp = -1d10`` to trace all markers. Codes of 101
and above are numerically unresolved orbits; they are excluded from the
confined fractions and reported in ``unresolved_fraction.dat``.

Step 5: Statistical error
-------------------------

Each marker is lost or confined independently, so the number of losses
:math:`k` out of :math:`N` resolved markers is binomial. The estimate
:math:`\hat f = k/N` has standard error

.. math::

   \sigma_f = \sqrt{\hat f (1-\hat f)/N}.

This plug-in error vanishes at :math:`k = 0` and then says nothing. For zero
observed losses the one-sided 95 % upper bound is

.. math::

   f < 1 - 0.05^{1/N} \simeq 3/N ,

which for :math:`N = 256` is 1.2 %. A loss fraction quoted as zero means "below
:math:`3/N`", not zero. For small :math:`k` use the Wilson interval
(:doc:`/user-guide/postprocessing`).

Save as ``tutorial1/plot_confined.py`` and run it in ``tutorial1``:

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np

   for case in ("qa", "qh"):
       t, passing, trapped, resolved = np.loadtxt(
           f"{case}/confined_fraction.dat", unpack=True)
       loss = 1.0 - (passing + trapped)
       n = resolved[-1]
       error = np.sqrt(loss[-1] * (1.0 - loss[-1]) / n)
       bound = 1.0 - 0.05 ** (1.0 / n)
       print(f"{case}: loss {loss[-1]:.3f} +- {error:.3f}"
             f" (zero-loss 95% bound {bound:.3f})")
       plt.semilogx(t, loss, label=case.upper())

   plt.xlabel("time (s)")
   plt.ylabel("loss fraction")
   plt.legend()
   plt.savefig("loss_fraction.png", dpi=150)

::

   qa: loss 0.121 +- 0.020 (zero-loss 95% bound 0.012)
   qh: loss 0.000 +- 0.000 (zero-loss 95% bound 0.012)

The QA/QH difference (12.1 ± 2.0 % against < 1.2 %) is significant. The plot
shows when the losses happen.

Step 6: Loss times of individual markers
----------------------------------------

``times_lost.dat`` has one line per marker (columns in
:doc:`/user-guide/outputs`): index, loss time in s, trapping parameter
:math:`\kappa`, start :math:`s`, ... Confined markers have loss time
``trace_time``, skipped markers ``-1``. The five earliest losses in QA:

.. code-block:: bash

   awk '$2 > 0 && $2 < 0.1' qa/times_lost.dat | LC_ALL=C sort -g -k2 | head -5

:math:`\kappa` (third column) is negative for passing and in :math:`(0,1)` for
trapped markers (:doc:`/model/guiding-centre`). In this QA run 28 of the 31
lost markers are trapped, and the first losses occur after about 20 µs.

What next
---------

* Vary ``sbeg``, ``ntestpart`` and ``trace_time`` and watch the loss fraction
  and its error.
* Add collisions with ``swcoll = .true.`` and background profiles
  (:doc:`/user-guide/running`).
* :doc:`wall-heat-loads` turns the QA ``results.nc`` into a heat-load pattern.
