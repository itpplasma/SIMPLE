SIMPLE Input Reference
======================

Generated from the Fortran declarations; do not edit this page by hand.
``simple.in`` requires ``&config`` and optionally contains ``&profiles``.
``&poincare`` belongs to the separate ``spectre_poincare.in`` input.
The readers are ``src/params.f90:read_config_namelist``,
``src/profiles.f90:read_profiles`` and ``app/spectre_poincare.f90:read_input``.

Validate an input from the repository root (Python standard library only)::

   python tools/check_input.py simple.in examples/*.in

The preflight reports unknown keys with suggestions, type errors, array
bounds and selected basic value constraints from the JSON schema. It returns
a nonzero status on errors. Runtime checks of files, fields and interacting
options remain the responsibility of SIMPLE. This is not a full Fortran parser.
It supports scalar and rank-one indexed/section assignments, repetition,
null values, comments, quoted strings, and slash or END group terminators.

Regenerate after configuring SIMPLE so that libneo sources are available::

   make configure CONFIG=Fast
   python tools/generate_input_reference.py
   python tools/generate_input_reference.py --check

For a custom libneo checkout, pass ``--libneo-source PATH`` to the generator.
CTest checks both this page and the :download:`JSON schema <input-reference.json>` for staleness.
Type selectors, shapes and initial defaults retain their source expressions;
named constants are also resolved in JSON's ``default_value`` and ``extent``.
A scalar initializer for an array applies to every element. Character lengths
are shown with the type. Defaults precede runtime alias resolution and depend
on the libneo sources selected for generation. An empty documentation cell
means no adjacent declaration comment was found; see ``examples/simple_full.in``
and ``DOC/config.md`` for additional usage, units and option interactions.

&config
-------

Source: ``src/params.f90:213``.

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 10 20 35

   * - Name
     - Type
     - Kind
     - Shape
     - Initial default
     - Declaration documentation
   * - ``notrace_passing``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``0``
     - Further configuration parameters
   * - ``nper``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``1000``
     -
   * - ``npoiper``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``100``
     -
   * - ``ntimstep``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``10000``
     -
   * - ``ntestpart``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``1024``
     -
   * - ``trace_time``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1d-1``
     -
   * - ``num_surf``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``1``
     -
   * - ``sbeg``
     - ``real``
     - ``dp``
     - ``1024``
     - ``0.5d0``
     -
   * - ``phibeg``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0d0``
     -
   * - ``thetabeg``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.0d0``
     -
   * - ``contr_pp``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``-1d0``
     -
   * - ``face_al``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1d0``
     -
   * - ``npoiper2``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``256``
     -
   * - ``n_e``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``2``
     -
   * - ``n_d``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``4``
     -
   * - ``netcdffile``
     - ``character(len=1000)``
     - ``default``
     - ``scalar``
     - ``''``
     - Empty = unset, apply 'wout.nc' default at usage
   * - ``ns_s``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``5``
     - <- spline order for R,Z,lambda over s
   * - ``ns_tp``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``5``
     - <- spline order for R,Z,lambda over theta and varphi
   * - ``multharm``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``5``
     -
   * - ``isw_field_type``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``0``
     -
   * - ``generate_start_only``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``startmode``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``1``
     -
   * - ``fusion_grid``
     - ``integer``
     - ``default``
     - ``3``
     - ``(/128, 32, 32/)``
     -
   * - ``grid_density``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0d0``
     -
   * - ``special_ants_file``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``integmode``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``EXPL_IMPL_EULER``
     -
   * - ``orbit_model``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``ORBIT_GC``
     -
   * - ``orbit_coord``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``0``
     - Chart for the full-orbit field+geometry: full orbit currently supports only orbit_coord = 1 (Boozer/chartmap), which shares the production GC field.
   * - ``relerr``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1d-13``
     -
   * - ``gpu_mode``
     - ``character(len=16)``
     - ``default``
     - ``scalar``
     - ``'off'``
     - GPU execution is configured in the same namelist as the CPU run. Native CUDA with tuned DOPRI is the default production path. The legacy OpenACC path remains available only when selected explicitly.
   * - ``gpu_backend``
     - ``character(len=16)``
     - ``default``
     - ``scalar``
     - ``'cuda-native'``
     -
   * - ``gpu_method``
     - ``character(len=16)``
     - ``default``
     - ``scalar``
     - ``'dopri'``
     -
   * - ``gpu_start_coordinates``
     - ``character(len=16)``
     - ``default``
     - ``scalar``
     - ``'reference'``
     -
   * - ``gpu_landreman``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.true.``
     -
   * - ``gpu_compact_init``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.false.``
     -
   * - ``gpu_work_order``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.true.``
     -
   * - ``gpu_num_devices``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``1``
     -
   * - ``gpu_t_block``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d-3``
     -
   * - ``gpu_loss_tau``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.1d0``
     -
   * - ``gpu_maxloss``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.02d0``
     -
   * - ``gpu_min_timestep``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d-10``
     -
   * - ``gpu_pilot_fraction``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.04d0``
     -
   * - ``gpu_particle_profile``
     - ``character(len=1000)``
     - ``default``
     - ``scalar``
     - ``''``
     -
   * - ``symplectic_newton_warning_mode``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.true.``
     -
   * - ``boundary_event_fraction_tolerance``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``-1d0``
     -
   * - ``boundary_event_radial_tolerance``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``-1d0``
     -
   * - ``canonical_grid_nr``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``62``
     -
   * - ``canonical_grid_ntheta``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``63``
     -
   * - ``canonical_grid_nphi``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``64``
     -
   * - ``canonical_ode_relerr``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1d-11``
     -
   * - ``tcut``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``-1d0``
     -
   * - ``nturns``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``8``
     -
   * - ``debug``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``class_plot``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     - <=AAA
   * - ``cut_in_per``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0d0``
     - <=AAA
   * - ``fast_class``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``vmec_b_scale``
     - ``real``
     - ``kind(1d0)``
     - ``scalar``
     - ``1d0``
     -
   * - ``vmec_rz_scale``
     - ``real``
     - ``kind(1d0)``
     - ``scalar``
     - ``1d0``
     -
   * - ``swcoll``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     - Colliding with D-T reactor plasma. TODO: Make configurable
   * - ``deterministic``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``old_axis_healing``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.True.``
     - DEPRECATED boolean switches (use axis_healing). Kept for back-compat and mapped to a mode with a deprecation warning; the next SIMPLE release will remove them and default axis_healing to 'polyfit'.
   * - ``old_axis_healing_boundary``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.True.``
     -
   * - ``axis_healing_power_law``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``rho_axis_heal``
     - ``real``
     - ``kind(1d0)``
     - ``scalar``
     - ``-1.0d0``
     - deprecated; use s_axis_heal
   * - ``axis_healing``
     - ``character(len=16)``
     - ``default``
     - ``scalar``
     - ``''``
     - Near-axis healing of the VMEC harmonics. Preferred selector: axis_healing = 'legacy'          - old default: discard min(\|m\|,4) surfaces and extrapolate. 'legacy_adaptive' - old adaptive boundary from determine_nheal_for_axis. 'powerlaw'        - rho\*\*\|m\| continuation below s_axis_heal. 'polyfit'         - rho\*\*\|m\| \* least-squares poly(s) below s_axis_heal. Recommended. The empty string means "unset": deprecated logical switches below are mapped to a mode with a warning (see resolve_axis_healing_mode).
   * - ``s_axis_heal``
     - ``real``
     - ``kind(1d0)``
     - ``scalar``
     - ``1.0d-2``
     - anchor flux for powerlaw/polyfit
   * - ``axis_healing_polyfit_degree``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``3``
     - polyfit polynomial degree in s
   * - ``am1``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``2.0d0``
     -
   * - ``am2``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``3.0d0``
     -
   * - ``z1``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d0``
     -
   * - ``z2``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d0``
     -
   * - ``densi1``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.5d14``
     -
   * - ``densi2``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.5d14``
     -
   * - ``tempi1``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d4``
     -
   * - ``tempi2``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d4``
     -
   * - ``tempe``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d4``
     -
   * - ``batch_size``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``2000000000``
     - Initialize large so batch mode is not default
   * - ``ran_seed``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``42``
     -
   * - ``reuse_batch``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``field_input``
     - ``character(len=1000)``
     - ``default``
     - ``scalar``
     - ``''``
     -
   * - ``coord_input``
     - ``character(len=1000)``
     - ``default``
     - ``scalar``
     - ``''``
     -
   * - ``wall_input``
     - ``character(len=1000)``
     - ``default``
     - ``scalar``
     - ``''``
     -
   * - ``wall_units``
     - ``character(len=16)``
     - ``default``
     - ``scalar``
     - ``'m'``
     -
   * - ``chart_boundary_kind``
     - ``character(len=16)``
     - ``default``
     - ``scalar``
     - ``'auto'``
     -
   * - ``integ_coords``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``-1000``
     - Sentinel: -1000 means user did not set it
   * - ``crossing_level``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``3``
     - SPECTRE interface treatment: only full orbit (3) is supported.
   * - ``spectre_fo_entry_larmor``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.5_dp``
     - Full orbit near SPECTRE interfaces: entry distance in radial gyroradii and the relative tolerance of RK45 on the sharp per-volume field. 1.5 keeps the switch's p_zeta invariant as at 3 (entry 1 breaks it) at about half the FO time (concept/spectre_loss_study/results/fo-residence).
   * - ``spectre_fo_tolerance``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0e-9_dp``
     -
   * - ``spectre_fo_grid``
     - ``integer``
     - ``default``
     - ``3``
     - ``(/32, 64, 32/)``
     - Per-volume FO spline points: rho, theta, zeta per field period.
   * - ``spectre_sbeg_is_toroidal_flux``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.false.``
     -
   * - ``spectre_ncon_r``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``48``
     - SPECTRE per-volume Meiss construction grid (n_r, n_th, n_phi). Each volume allocates rank-3 arrays plus batch splines of this size, so the total peak memory scales as Mvol\*n_r\*n_th\*n_phi. spectre_ncon_phi = -1 means automatic: 32 for fields with toroidal harmonics and a minimal phi grid for axisymmetric fields (all n = 0), the dominant memory saver for tokamak cases. A positive value forces that phi count verbatim, bypassing the axisymmetric clamp.
   * - ``spectre_ncon_th``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``48``
     -
   * - ``spectre_ncon_phi``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``-1``
     -
   * - ``spectre_ncon_order``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``5``
     -
   * - ``spectre_ncon_ode_max_steps``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``1000000``
     -
   * - ``spectre_ncon_ode_relerr``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0e-2_dp``
     -
   * - ``output_results_netcdf``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``output_error``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``output_orbits_macrostep``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``macrostep_time_grid``
     - ``character(len=16)``
     - ``default``
     - ``scalar``
     - ``'linear'``
     -
   * - ``checkpoint_interval``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``10.0d0``
     - Legacy progress/checkpoint interval, retained for input compatibility. Tracing performs no timed I/O; results and event totals are emitted after the OpenMP team joins.
   * - ``restart``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.False.``
     -
   * - ``local``
     - ``logical``
     - ``default``
     - ``scalar``
     - ``.true.``
     - Removed &config keys. They are still read so that old input files run; a value other than the historical default prints a warning.
   * - ``loopskip``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``0``
     -
   * - ``bmod_ref``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``5.0d4``
     -
   * - ``macrostep_log_ratio``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d3``
     -

&profiles
---------

Source: ``src/profiles.f90:36``.

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 10 20 35

   * - Name
     - Type
     - Kind
     - Shape
     - Initial default
     - Declaration documentation
   * - ``profile_type``
     - ``character(len=32)``
     - ``default``
     - ``scalar``
     - ``"power_series"``
     -
   * - ``te_scale``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d4``
     -
   * - ``ti1_scale``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d4``
     -
   * - ``ti2_scale``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d4``
     -
   * - ``ni1_scale``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.5d20``
     -
   * - ``ni2_scale``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.5d20``
     -
   * - ``te_coef``
     - ``real``
     - ``dp``
     - ``MAX_POWER_SERIES``
     - ``[1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]``
     -
   * - ``ti1_coef``
     - ``real``
     - ``dp``
     - ``MAX_POWER_SERIES``
     - ``[1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]``
     -
   * - ``ti2_coef``
     - ``real``
     - ``dp``
     - ``MAX_POWER_SERIES``
     - ``[1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]``
     -
   * - ``ni1_coef``
     - ``real``
     - ``dp``
     - ``MAX_POWER_SERIES``
     - ``[1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]``
     -
   * - ``ni2_coef``
     - ``real``
     - ``dp``
     - ``MAX_POWER_SERIES``
     - ``[1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]``
     -
   * - ``te_p1``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d0``
     -
   * - ``te_p2``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.0d0``
     -
   * - ``ti1_p1``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d0``
     -
   * - ``ti1_p2``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.0d0``
     -
   * - ``ti2_p1``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d0``
     -
   * - ``ti2_p2``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.0d0``
     -
   * - ``ni1_p1``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d0``
     -
   * - ``ni1_p2``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.0d0``
     -
   * - ``ni2_p1``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``1.0d0``
     -
   * - ``ni2_p2``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.0d0``
     -

&poincare
---------

Source: ``app/spectre_poincare.f90:32``.

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 10 20 35

   * - Name
     - Type
     - Kind
     - Shape
     - Initial default
     - Declaration documentation
   * - ``spectre_file``
     - ``character(len=1024)``
     - ``default``
     - ``scalar``
     - ``'spectre.h5'``
     -
   * - ``zeta0``
     - ``real``
     - ``dp``
     - ``scalar``
     - ``0.0_dp``
     -
   * - ``nturns``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``200``
     -
   * - ``nsteps_per_period``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``256``
     -
   * - ``nseeds_per_volume``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``4``
     -
   * - ``volumes``
     - ``integer``
     - ``default``
     - ``max_volumes``
     - ``-1``
     -
   * - ``seed_file``
     - ``character(len=1024)``
     - ``default``
     - ``scalar``
     - ``''``
     -
   * - ``spline_n_s``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``32``
     -
   * - ``spline_n_theta``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``32``
     -
   * - ``spline_n_zeta``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``32``
     -
   * - ``spline_order``
     - ``integer``
     - ``default``
     - ``scalar``
     - ``3``
     -
