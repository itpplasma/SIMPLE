Outputs
=======

This page lists every file SIMPLE writes, its columns and the meaning of the
exit codes. Times are in s; angles in radians, not wrapped; positions in the
reference coordinates (:math:`s` for VMEC, :math:`\rho = \sqrt{s}` for Boozer
chartmaps, the native radius for other chartmaps).

Phase-space state
-----------------

=====  ==================================  ============================================
Index  Symbol                              Meaning
=====  ==================================  ============================================
1      :math:`s` or :math:`\rho`           radial reference coordinate
2      :math:`\vartheta`                   poloidal reference angle
3      :math:`\varphi`                     toroidal reference angle
4      :math:`v/v_0`                       speed over birth speed; 1 without collisions
5      :math:`\lambda = v_\parallel/v`     pitch
=====  ==================================  ============================================

Text files
----------

``confined_fraction.dat``
   One row per macrostep: time, confined fraction of passing markers,
   confined fraction of trapped markers, number of resolved markers. Both
   fractions are divided by the resolved markers at that time; the loss
   fraction is :math:`1 - f_\text{pass} - f_\text{trap}`. Passing means
   :math:`\kappa < 0` at the start.
``unresolved_fraction.dat``
   One row per macrostep: time, fraction of unresolved markers, ``ntestpart``.
``times_lost.dat``
   One row per marker: index, loss time, trapping parameter :math:`\kappa`,
   start radius, perpendicular invariant :math:`J_\perp`, end state (5
   columns). The loss time is ``trace_time`` for confined markers, ``-1``
   for markers not traced (exit code 3) and ``NaN`` for unresolved ones.
``orbit_exit_code.dat``
   One row per marker: index, exit code, loss time, radial residual and time
   width of the located boundary crossing.
``start.dat``
   Start states, one per row (written by ``startmode`` 1, 5, 7; read by 2, 4,
   6).
``avg_inverse_t_lost.dat``
   Mean of :math:`1/t_\text{loss}` over lost markers (0 if none).
``class_parts.dat``
   With classification (``tcut > 0`` or ``class_plot``): index, start
   radius, :math:`J_\perp`, three class codes (:doc:`classification`).
``bminmax.dat``
   With classification and ``num_surf > 1``: :math:`s`,
   :math:`B_\text{min}(s)`, :math:`B_\text{max}(s)` in G.
``alpha_source.dat``
   Source normalisation of ``startmode = 7`` (:doc:`/alpha-birth`).

:math:`\kappa` and :math:`J_\perp` are defined in
:doc:`/model/guiding-centre`.

Exit codes
----------

======  ======================================================================
Code    Meaning
======  ======================================================================
0       traced to ``trace_time``: confined
1       lost through the LCFS (or a chartmap boundary declared ``'lcfs'``)
2       lost to the wall (STL wall or chartmap boundary declared ``'wall'``)
3       not traced (``contr_pp``, ``notrace_passing``); counted confined
4       all recovery methods failed for a collisionless marker in the core
        (:math:`s < 0.01`); counted confined
101     unresolved: left the coordinate domain (``chart_boundary_kind = 'domain'``)
102     unresolved: Newton iteration did not converge
103     unresolved: linear solve failed
104     unresolved: boundary crossing could not be located
105     unresolved: full-orbit step failed
======  ======================================================================

Unresolved markers (101–105) are excluded from the confined fractions and
counted in ``unresolved_fraction.dat``. A large unresolved fraction means a
too large step (``npoiper2``), too coarse field splines or a bad equilibrium
near the axis; it is not a physical loss.

NetCDF results
--------------

``output_results_netcdf = .true.`` writes ``results.nc`` with dimension
``particle`` (``ntestpart``) and global attributes holding the run settings
(``trace_time``, ``ntestpart``, ``integmode``, ``cart_units``, input files,
...).

======================================  ===============================================
Variable                                Meaning
======================================  ===============================================
``zstart``, ``zend``                    start and end state (``particle``, ``phase``);
                                        shape (n, 5) in Python
``xstart_cart``, ``xend_cart``          Cartesian start and end position (``units``)
``times_lost``, ``orbit_exit_code``     as in the text files
``boundary_event_radial_residual``      column 4 of ``orbit_exit_code.dat``
``boundary_event_time_width``           column 5 of ``orbit_exit_code.dat``
``trap_par``, ``perp_inv``              :math:`\kappa`, :math:`J_\perp`
``iclass``                              three class codes; zero without classification
``class_lost``                          loss flag of the classification run
``wall_hit``, ``wall_hit_cart``         STL wall hit flag and intersection point
``wall_hit_triangle``                   hit STL facet, 0-based; only with ``wall_input``
``wall_hit_normal_cart``                facet normal at the hit
``wall_hit_cos_incidence``              incidence cosine; 1 is normal incidence
``wall_hit_angle_rad``                  incidence angle; 0 is normal incidence
======================================  ===============================================

``output_orbits_macrostep = .true.`` writes ``orbits.nc`` with every
marker's state at every macrostep: ``time``, ``s``, ``theta``, ``phi``,
``p_abs`` (:math:`v/v_0`), ``v_par`` (:math:`\lambda`) and, when a reference
geometry exists, ``R`` and ``Z``. The attribute ``radial_coordinate`` names
the meaning of ``s``.

All text files load with ``numpy.loadtxt``. Plots: :doc:`postprocessing`;
``pysimple.plotting.plot_loss_statistics`` adds energy-resolved loss panels.
