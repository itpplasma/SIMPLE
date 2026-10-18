Running SIMPLE
==============

Command line
------------

.. code-block:: bash

   simple.x [config_file]

``simple.x`` reads ``simple.in`` from the working directory, or the file given
as the first argument, and writes all results to the working directory. The
equilibrium and other input files are found relative to the working directory
unless their paths are absolute. The Python interface runs the same code
(:doc:`python`).

Parallelisation
---------------

Particles are independent and are distributed over OpenMP threads. Set the
number of threads with ``OMP_NUM_THREADS``; on machines with several sockets,
``OMP_PROC_BIND=spread`` and ``OMP_PLACES=cores`` keep the threads on separate
cores. The field setup before the tracing is partly serial and takes a few
seconds to a minute, depending on the field source and resolution. The log
reports the time of each phase. GPU runs are described in :doc:`gpu`, and
cluster job scripts for Perlmutter, aCluster and sCluster are in ``RUN/``.

Time and particle settings
--------------------------

``trace_time``
   Physical tracing time in seconds. Particles that are not lost before this
   time count as confined.
``ntestpart``
   Number of particles (markers). The statistical error of a loss fraction
   :math:`f` is :math:`\sqrt{f(1-f)/N}`.
``ntimstep``
   Number of output times (macrosteps) between 0 and ``trace_time``; one line
   of ``confined_fraction.dat`` per macrostep. ``macrostep_time_grid = 'log'``
   spaces them logarithmically, which resolves the early losses better.
``npoiper2``
   Number of integration steps (microsteps) per path length
   :math:`2\pi R_0` at the birth speed. The default 256 is safe for
   symplectic Euler; see ``DOC/config.md`` in :doc:`/reference/design-notes`
   for when fewer steps are acceptable.
``contr_pp``, ``notrace_passing``
   Passing particles with trapping parameter below ``contr_pp`` (default
   :math:`-1`) are counted as confined without tracing, because their orbits
   stay close to a flux surface. ``contr_pp = -1d10`` traces all particles.
   ``notrace_passing = 1`` skips all passing particles.
``deterministic``, ``ran_seed``
   A fixed seed for sampling and collisions, for reproducible runs.

Particle species and energy
---------------------------

The default particle is a 3.5 MeV alpha particle. ``n_e`` and ``n_d`` set the
charge and mass number, and ``facE_al`` divides the energy:
:math:`E = 3.5\,\text{MeV}/\text{facE\_al}`. For example, 100 keV protons are
``n_e = 1``, ``n_d = 1d0``, ``facE_al = 35d0``.

Integrators
-----------

``integmode`` selects the integrator:

.. list-table::
   :header-rows: 1
   :widths: 10 35 15 40

   * - Value
     - Method
     - Order
     - Notes
   * - 1
     - Symplectic Euler (explicit-implicit)
     - 1
     - Default; fastest; first order but symplectic.
   * - 2
     - Symplectic Euler (implicit-explicit)
     - 1
     - Adjoint of 1.
   * - 3
     - Implicit midpoint
     - 2
     - Symmetric and symplectic.
   * - 4--7
     - Gauss-Legendre collocation, 1 to 4 stages
     - 2, 4, 6, 8
     - High accuracy per step; one Newton solve for all stages.
   * - 15
     - Lobatto IIIA-IIIB, 3 stages
     - 4
     - Partitioned symplectic Runge-Kutta.
   * - 0
     - Adaptive Runge-Kutta 4/5
     - adaptive
     - Not symplectic; works with all field types.

The symplectic integrators (``integmode > 0``) need canonical coordinates,
``integ_coords`` 0, 2, 3, 4 or 6. Their implicit equations are solved by
Newton iterations to the relative tolerance ``relerr``. The methods are
described in :doc:`/model/integrators`.

Collisions
----------

``swcoll = .true.`` adds Coulomb collisions of the test particles with a
background of electrons and two ion species. The background is defined by
radial profiles in the ``&profiles`` namelist, as power series or as
:math:`f(s) = f_0 (1 - s^{p_1})^{p_2}`:

.. code-block:: fortran

   &config
   swcoll = .true.
   am1 = 2d0, am2 = 3d0          ! D and T
   Z1 = 1d0, Z2 = 1d0
   /
   &profiles
   profile_type = 'two_power'
   Te_scale = 1.5d4, Te_p1 = 2d0, Te_p2 = 1d0         ! eV
   Ti1_scale = 1.5d4, Ti1_p1 = 2d0, Ti1_p2 = 1d0
   Ti2_scale = 1.5d4, Ti2_p1 = 2d0, Ti2_p2 = 1d0
   ni1_scale = 5d19, ni1_p1 = 2d0, ni1_p2 = 0.5d0     ! m^-3
   ni2_scale = 5d19, ni2_p1 = 2d0, ni2_p2 = 0.5d0
   /

With collisions, ``trace_time`` should cover the slowing-down time, and the
energy carried to the wall is reduced (column 4 of the final state,
:math:`v/v_0`). Collisions cannot be combined with orbit classification or
full-orbit tracing. The model is described in :doc:`/model/guiding-centre`.

Restart and batches
-------------------

``restart = .true.`` continues a run from ``times_lost.dat`` and
``orbit_exit_code.dat``; particles that reached ``trace_time`` are not traced
again. ``batch_size`` smaller than ``ntestpart`` traces a random subset of
the start positions, whose indices are written to ``batch.dat``.

Numerical robustness
--------------------

``symplectic_newton_warning_mode = .true.`` (default) accepts bounded Newton
corrections when the iteration limit is reached and retries a failed step with
smaller symplectic steps and then with the adaptive Runge-Kutta method. Orbits
that fail anyway end with a numerical exit code (101 to 105) and a loss time
``NaN``; they are reported as unresolved and not as lost
(:doc:`outputs`). ``DOC/config.md`` (:doc:`/reference/design-notes`) gives
the details and the related tolerances.
