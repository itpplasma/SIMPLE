Quick start
===========

This page computes the alpha-particle loss fraction of a public
quasi-axisymmetric reactor-scale VMEC equilibrium (Landreman and Paul 2022,
from the simsopt test files) with ``examples/simple.in``.

Run
---

From the repository root, after :doc:`installation`:

.. code-block:: bash

   mkdir -p run
   cp examples/simple.in run/
   curl -L -o run/wout.nc https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc
   cd run
   ../build/simple.x simple.in

.. literalinclude:: ../../examples/simple.in
   :language: fortran
   :caption: examples/simple.in

The run traces 128 alpha particles of 3.5 MeV for 10 ms, started on the
surface :math:`s = 0.3` (normalised toroidal flux) with random angles and
uniformly random pitch :math:`\lambda = v_\parallel/v`. ``isw_field_type = 2``
integrates in Boozer coordinates. All other keys keep the defaults of
``examples/simple_full.in``. The run takes 20–30 s on a workstation; most of
it is field setup, which the log reports phase by phase.

Result
------

.. code-block:: bash

   tail -1 confined_fraction.dat

::

   1.0000000000000000E-002  0.70312500000000000  0.29687500000000000  128

The columns are time (s), confined fraction of passing markers, confined
fraction of trapped markers and the number of numerically resolved markers.
The loss fraction at 10 ms is :math:`1 - f_\text{pass} - f_\text{trap} = 0`:
no marker left the plasma. Starts are random, so another run differs; in this
optimised configuration 0 to a few of 128 markers are lost. Add
``deterministic = .true.`` for a fixed seed.

``orbit_exit_code.dat`` shows code 3 for about half of the markers: deeply
passing markers (trapping parameter :math:`\le` ``contr_pp``) are counted as
confined without tracing. ``times_lost.dat`` holds the loss time and end state
of every marker. :doc:`/user-guide/outputs` describes all files and codes.

Next
----

* :doc:`/tutorials/loss-fraction`: QA versus QH with statistical errors.
* :doc:`/user-guide/running`: markers, time, integrator, starts, collisions.
* :doc:`/model/index`: what is computed.
