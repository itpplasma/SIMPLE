Quick start
===========

This page runs a first loss computation in a few minutes. It uses a public
quasi-axisymmetric reactor-scale VMEC equilibrium (Landreman and Paul 2022)
from the simsopt test files and the input ``examples/simple.in``.

Run
---

From the repository root, after :doc:`installation`:

.. code-block:: bash

   mkdir -p run
   cp examples/simple.in run/simple.in
   wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc \
       -O run/wout.nc
   cd run
   ../build/simple.x simple.in

The input file is a Fortran namelist:

.. literalinclude:: ../../examples/simple.in
   :language: fortran
   :caption: examples/simple.in

It traces 128 alpha particles with 3.5 MeV for 10 ms. They start on the flux
surface with normalised toroidal flux :math:`s = 0.3` with random pitch.
``isw_field_type = 2`` selects Boozer coordinates, in which SIMPLE builds the
canonical guiding-centre equations. All other keys keep their defaults, which
are listed in ``examples/simple_full.in`` (:doc:`/reference/input-parameters`).

On a workstation the run takes about ten seconds. Most of it is the setup of
the field splines and the Boozer transformation. The log ends with lines such
as::

   INFO: Particle tracing completed       2.902 s (Total:       10.572 s)
    Total field evaluations:              74735691
   INFO: Output writing completed       0.026 s (Total:       10.598 s)

Results
-------

The run directory now contains the result files. The most important one is
``confined_fraction.dat``. Its columns are the time in seconds, the confined
fraction of passing particles, the confined fraction of trapped particles and
the number of numerically resolved particles:

.. code-block:: bash

   tail -1 confined_fraction.dat

::

   1.0000000000000000E-002  0.76562500000000000  0.22656250000000000  128

The loss fraction is :math:`1 - (f_\text{pass} + f_\text{trap})` at the last
time, here :math:`1 - 0.992 = 0.008`: one of 128 particles left the plasma
within 10 ms. The start positions and pitches are random, so another run
gives slightly different numbers; with 128 particles, zero to a few losses are
expected in this optimised configuration. Add ``deterministic = .true.`` to
``simple.in`` for a fixed random seed. ``times_lost.dat`` contains the loss
time of each particle, and :doc:`/user-guide/outputs` describes all files.

The script ``examples/run_example.sh`` runs the same case in
``/tmp/simple_example``.

Next steps
----------

* :doc:`/tutorials/loss-fraction` computes and compares the loss fractions of
  a quasi-axisymmetric and a quasi-helical equilibrium and plots the result.
* :doc:`/user-guide/running` explains the choice of particle number, tracing
  time, integrator and starting positions.
* :doc:`/model/index` explains what is computed.
