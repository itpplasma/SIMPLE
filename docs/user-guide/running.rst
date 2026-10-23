Running SIMPLE
==============

This page covers the settings of a run. Field inputs are in
:doc:`field-sources`, outputs in :doc:`outputs`.

Method: :doc:`/model/integrators`.

Command line and threads
------------------------

.. code-block:: bash

   OMP_NUM_THREADS=16 simple.x [config_file]

``simple.x`` reads ``simple.in`` or the given file, resolves relative input
paths against the working directory and writes all results there. Markers are
distributed over OpenMP threads. On multi-socket nodes,
``OMP_PROC_BIND=spread OMP_PLACES=cores`` keeps threads on separate cores.
Field setup before tracing is partly serial; the log reports the time of
each phase. GPU runs: :doc:`gpu`; cluster jobs: :doc:`/hpc/index`.

Time and markers
----------------

============================  ==============  ==============================================
Key                           Default         Meaning
============================  ==============  ==============================================
``trace_time``                ``1d-1``        tracing time in s; markers not lost by then
                                              count as confined
``ntestpart``                 1024            number of markers :math:`N`; a loss fraction
                                              :math:`f` has standard error
                                              :math:`\sqrt{f(1-f)/N}`
``ntimstep``                  10000           output times (macrosteps) up to
                                              ``trace_time``, rows of ``confined_fraction.dat``
``macrostep_time_grid``       ``'linear'``    ``'log'``: logarithmic spacing (last/first
                                              interval 10), resolves early losses
``npoiper2``                  256             steps per path length :math:`2\pi R_0` at
                                              speed :math:`v_0`
``contr_pp``                  ``-1d0``        markers with :math:`\kappa \le` ``contr_pp``
                                              (deeply passing) count as confined untraced
                                              (exit code 3); ``-1d10`` traces all
``notrace_passing``           0               1: all passing markers count as confined
``deterministic``             ``.false.``     fixed seed ``ran_seed`` (default 42) for
                                              sampling and collisions
============================  ==============  ==============================================

Confidence intervals: :doc:`postprocessing`. Fewer steps:
:doc:`/reference/numerical-settings`. :math:`\kappa` (trapping parameter):
:doc:`/model/guiding-centre`.

Species and energy
------------------

``n_e`` (charge number, default 2), ``n_d`` (mass number, default ``4d0``)
and ``facE_al`` (default ``1d0``, energy :math:`E_0 = 3.5\,\text{MeV}/`
``facE_al``) set the species: 100 keV protons are ``n_e = 1, n_d = 1d0,
facE_al = 35d0``. Alternatively set ``particle_energy_eV = 100000``
with the default ``facE_al = 1``. These species settings apply to both
orbit integration and full Coulomb collisions.

Start positions
---------------

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - ``startmode``
     - Start positions
   * - 1 (default)
     - Random points on one field line on the surface ``sbeg(1)`` (default
       0.5), started at ``thetabeg``, ``phibeg`` and followed for ``nper``
       field periods; volume-weighted. ``grid_density`` in (0, 1) places
       :math:`(\text{int}(1/\text{grid\_density} - 1))^2` markers on a
       :math:`(\vartheta, \varphi)` grid instead and overrides ``ntestpart``.
   * - 2
     - Read ``start.dat``: :math:`(s, \vartheta, \varphi, v/v_0, \lambda)`
       per line, reference coordinates.
   * - 3
     - Read an ANTS-format file (``start_ants.dat`` with
       ``special_ants_file = .true.``, else ``start.dat``).
   * - 4
     - Random batch of ``batch_size`` lines of ``start.dat``.
   * - 5
     - Uniform random :math:`s` in the range set by ``num_surf`` (0:
       :math:`[0,1]`, 1: :math:`[0,` ``sbeg(1)`` :math:`]`, 2:
       :math:`[` ``sbeg(1)``, ``sbeg(2)`` :math:`]`), clamped to
       :math:`s \ge 0.01`; random angles.
   * - 6
     - Read ``start.dat`` with positions in integration coordinates.
   * - 7
     - D-T fusion births from ``&profiles``: :doc:`/alpha-birth`.

Modes 1, 5 and 7 write ``start.dat``; ``generate_start_only = .true.`` stops
after that. Generated markers have :math:`v/v_0 = 1` and pitch
:math:`\lambda = v_\parallel/v` uniform in :math:`[-1, 1]`. Angles are in
radians and are not wrapped.

Integrators
-----------

=============  ======================================  ==========  ==================================
``integmode``  Method                                  Order       Notes
=============  ======================================  ==========  ==================================
1 (default)    symplectic Euler, explicit-implicit     1           cheapest per step
2              symplectic Euler, implicit-explicit     1           adjoint of 1
3              implicit midpoint                       2           symmetric
4–7            Gauss-Legendre, 1–4 stages              2, 4, 6, 8  one Newton solve for all stages
15             Lobatto IIIA-IIIB, 3 stages             4           partitioned
0              adaptive Runge-Kutta 4(5)               adaptive    not symplectic; any ``integ_coords``
-1             as 0                                    adaptive    skips loading ``field_input``
=============  ======================================  ==========  ==================================

``integmode > 0`` needs canonical coordinates (``integ_coords`` not 1 or 5).
Implicit steps are solved by Newton iteration to relative tolerance
``relerr`` (default ``1d-13``), which also sets the Runge-Kutta tolerance.

``symplectic_newton_warning_mode = .true.`` (default) accepts bounded Newton
corrections at the iteration limit and retries a failed step with smaller
symplectic steps, then with Runge-Kutta. Markers that still fail get a
numerical exit code (101–105), loss time ``NaN`` and count as unresolved, not
lost. ``.false.`` ends the marker at the first failure. Tolerances:
:doc:`/reference/numerical-settings`.

Collisions
----------

``swcoll = .true.`` adds Coulomb collisions with electrons and two ion
species (``am1``, ``am2``, ``Z1``, ``Z2``; default D and T). Density and
temperature profiles :math:`f(s)` come from ``&profiles``, either
``'power_series'`` :math:`f = f_0 \sum_{n=1}^{10} c_n s^{n-1}` or
``'two_power'`` :math:`f = f_0 (1 - s^{p_1})^{p_2}`:

.. code-block:: fortran

   &config
   swcoll = .true.
   am1 = 2d0, am2 = 3d0, Z1 = 1d0, Z2 = 1d0     ! D and T
   /
   &profiles
   profile_type = 'two_power'
   Te_scale  = 1.5d4, Te_p1  = 2d0, Te_p2  = 1d0    ! eV
   Ti1_scale = 1.5d4, Ti1_p1 = 2d0, Ti1_p2 = 1d0
   Ti2_scale = 1.5d4, Ti2_p1 = 2d0, Ti2_p2 = 1d0
   ni1_scale = 5d19,  ni1_p1 = 2d0, ni1_p2 = 0.5d0  ! m^-3
   ni2_scale = 5d19,  ni2_p1 = 2d0, ni2_p2 = 0.5d0
   /

Limits that change results:

* ``n_d`` and ``n_e`` select the test-ion mass and charge. Energy is set by
  ``particle_energy_eV`` or ``facE_al``.
* Collisions cannot be combined with classification (``tcut > 0``,
  ``class_plot``, ``fast_class``), full orbits or GPU tracing.
* ``trace_time`` must cover the slowing-down time for complete heat loads.
  The marker energy at loss is :math:`(v/v_0)^2 E_0`, from column 4 of the end
  state.

Collision operator: :doc:`/model/collisions`.

Restart and batches
-------------------

``restart = .true.`` reads ``times_lost.dat`` and ``orbit_exit_code.dat``;
markers with exit code 0, 1, 2 or 4 keep their result, the others are traced
again. ``batch_size`` smaller than ``ntestpart`` traces a random subset of
the start positions, whose indices go to ``batch.dat`` (``reuse_batch =
.true.`` reads them back).
