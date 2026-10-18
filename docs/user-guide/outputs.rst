Outputs and their meaning
=========================

SIMPLE writes text files in every run and NetCDF files on request. Times are
in seconds, angles in radians and not wrapped to :math:`[0, 2\pi)`, and
positions in the reference coordinates (:math:`s` for VMEC, :math:`\rho` for
chartmaps).

Phase-space state
-----------------

Start and end states have five components:

.. list-table::
   :header-rows: 1
   :widths: 10 25 65

   * - Index
     - Symbol
     - Meaning
   * - 1
     - :math:`s` or :math:`\rho`
     - Radial reference coordinate: normalised toroidal flux for VMEC,
       chartmap radius for chartmaps.
   * - 2
     - :math:`\vartheta`
     - Poloidal angle of the reference coordinates.
   * - 3
     - :math:`\varphi`
     - Toroidal angle of the reference coordinates.
   * - 4
     - :math:`v/v_0`
     - Speed normalised to the birth speed; 1 without collisions.
   * - 5
     - :math:`\lambda = v_\parallel/v`
     - Pitch.

Text files
----------

``confined_fraction.dat``
   One line per macrostep: time, confined fraction of passing particles,
   confined fraction of trapped particles, number of resolved particles. The
   fractions are normalised to the resolved particles at that time. The loss
   fraction is :math:`1 - f_\text{pass} - f_\text{trap}`.
``unresolved_fraction.dat``
   One line per macrostep: time, fraction of numerically unresolved orbits,
   total particle number.
``times_lost.dat``
   One line per particle: index, loss time, trapping parameter, start radius,
   perpendicular invariant :math:`J_\perp`, final state (5 columns). The loss
   time is ``trace_time`` for confined particles, ``-1`` for particles that
   were not traced (``contr_pp``) and ``NaN`` for unresolved orbits.
``orbit_exit_code.dat``
   One line per particle: index, exit code, loss time, radial residual and
   time width of the located boundary event.
``start.dat``
   Start states, one per line (written for ``startmode`` 1 and 5, read for 2
   and 4).
``avg_inverse_t_lost.dat``
   Mean of the inverse loss time of lost particles.
``class_parts.dat``, ``bminmax.dat``
   Classification results and the :math:`B_\text{min}`, :math:`B_\text{max}`
   table, when classification is on (:doc:`classification`).

Exit codes
----------

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Code
     - Meaning
   * - 0
     - Traced to ``trace_time`` (confined).
   * - 1
     - Lost: crossed the last closed flux surface.
   * - 2
     - Lost: hit the wall (STL wall or chartmap boundary declared as wall).
   * - 3
     - Not traced; counted as confined (``contr_pp``, ``notrace_passing``).
   * - 4
     - Numerically confined: all recovery methods failed in the core
       (:math:`s < 0.01`) of a collisionless orbit; counted as confined.
   * - 101
     - Numerical: left the domain of the coordinates (``chart_boundary_kind =
       'domain'``).
   * - 102
     - Numerical: Newton iteration did not converge.
   * - 103
     - Numerical: linear solve failed.
   * - 104
     - Numerical: boundary event could not be located.
   * - 105
     - Numerical: full-orbit step failed.

Codes 101 to 105 are unresolved orbits. They are excluded from the confined
fractions and reported in ``unresolved_fraction.dat``. A large unresolved
fraction indicates a too large time step, a too coarse field resolution or a
problem of the equilibrium near the axis.

Trapping parameter
------------------

The trapping parameter of a particle is evaluated at its start point,

.. math::

   \kappa = \left(\frac{\mu B_\text{max}}{E} - 1\right)
            \frac{B_\text{min}}{B_\text{max} - B_\text{min}}
          = \left((1 - \lambda^2)\frac{B_\text{max}}{B} - 1\right)
            \frac{B_\text{min}}{B_\text{max} - B_\text{min}},

with :math:`B_\text{min}`, :math:`B_\text{max}` on the start surface.
Negative values are passing, 0 is the trapped-passing boundary and 1 is a
deeply trapped particle reflecting at :math:`B_\text{min}`. The split of
``confined_fraction.dat`` into passing and trapped uses :math:`\kappa < 0`.

NetCDF results
--------------

``output_results_netcdf = .true.`` writes ``results.nc`` with one entry per
particle (dimension ``particle``):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable
     - Meaning
   * - ``zstart``, ``zend``
     - Start and end state, shape (``particle``, ``phase``).
   * - ``xstart_cart``, ``xend_cart``
     - Cartesian start and end position (unit in the ``units`` attribute).
   * - ``times_lost``, ``orbit_exit_code``
     - As in the text files.
   * - ``trap_par``, ``perp_inv``
     - Trapping parameter, perpendicular invariant.
   * - ``iclass``
     - Classification results (3 per particle).
   * - ``class_lost``
     - Loss flag of the classification run; it is only set when
       classification is on.
   * - ``wall_hit``, ``wall_hit_cart``
     - STL wall hit flag and intersection point.
   * - ``wall_hit_normal_cart``, ``wall_hit_cos_incidence``,
       ``wall_hit_angle_rad``
     - Wall normal and incidence angle at the hit (0 is normal incidence).

Global attributes record the run settings, for example ``trace_time`` and
``ntestpart``. In Python, the array ``zend`` has the shape (particles, 5);
transpose it to index the phase-space component first.

``output_orbits_macrostep = .true.`` writes ``orbits.nc`` with the state of
every particle at every macrostep (variables ``s``, ``theta``, ``phi``,
``p_abs``, ``v_par`` and, when a reference geometry exists, ``R`` and ``Z``).
The attribute ``radial_coordinate`` gives the meaning of ``s``.

Reading results in Python
-------------------------

``pysimple.plotting`` reads a run directory and plots the standard loss
statistics:

.. code-block:: python

   from pysimple.plotting import load_loss_data, plot_loss_statistics

   data = load_loss_data("run")
   plot_loss_statistics(data, "loss_statistics.png", show=False)

``examples/plot_losses.py`` and ``examples/plot_energy_loss.py`` are
command-line versions.
