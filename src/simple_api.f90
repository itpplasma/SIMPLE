module simple_api
    !> Library interface to one SIMPLE run.
    !>
    !> The run configuration is the &config namelist held in module `params`;
    !> set it with simple_read_config and/or by assigning `params` variables,
    !> then call simple_init. The remaining calls repeat the steps of simple.x
    !> (sampling, tracing, classification, output) and add single-orbit
    !> diagnostics. All states are (s, theta, phi, v/v0, v_par/v) in reference
    !> coordinates. The Python package pysimple wraps this module and `params`.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use params, only: ntestpart, startmode, num_surf, sbeg, zstart, zend, &
        times_lost, trap_par, perp_inv, iclass, class_passing, class_lost, &
        ntimstep, bmod00, bmin, bmax, reallocate_arrays, reset_config, &
        finalize_config, isw_field_type, integ_coords, field_input, config
    use simple_profiles, only: read_profiles, reset_profiles, set_profile_table
    use magfie_sub, only: TEST, CANFLUX, BOOZER
    use field_boozer_chartmap, only: is_boozer_chartmap
    use field_can_mod, only: integ_to_ref
    use simple_main, only: validate_orbit_model_config, init_run, &
        init_sampling_field, sample_start, init_tracing_field, init_counters, &
        trace_parallel, classify_parallel, trace_orbit, write_output, &
        compute_canonical_frequencies_flat, invariants_from_state_flat, &
        potato_invariants_to_simple_flat, states_from_invariants_flat, &
        trace_to_cut_flat

    implicit none
    private

    public :: simple_reset_config, simple_read_config, simple_set_config, &
              simple_init, &
              simple_sample, simple_sample_surface, simple_sample_volume, &
              simple_set_start, simple_trace, simple_classify, simple_get_results, &
              simple_get_classification, simple_write_output, simple_trace_orbit, &
              simple_orbit_frequencies, simple_invariants_from_state, &
              simple_states_from_invariants, simple_potato_invariants, &
              simple_trace_to_cut, simple_reference_info, simple_wall_supported, &
              simple_cylindrical_to_reference, simple_reference_to_cylindrical, &
              simple_set_profile_table

    integer, parameter :: FIELD_NONE = 0, FIELD_SAMPLING = 1, FIELD_TRACING = 2

    type(tracer_t), save :: tracer
    integer, save :: active_field = FIELD_NONE
    real(dp), save :: sampling_surface = -1.0_dp
    !> reallocate_arrays resets the starting-line arrays (xstart, volstart) that
    !> init_starting_surf fills, so they must be integrated again before sampling.
    logical, save :: starting_line_stale = .true.

contains

    logical function simple_wall_supported()
        use stl_wall_intersection, only: stl_wall_available

        simple_wall_supported = stl_wall_available
    end function simple_wall_supported

    subroutine simple_set_profile_table(n, table, status)
        integer, intent(in) :: n
        real(dp), intent(in) :: table(6, n)
        integer, intent(out) :: status

        call set_profile_table(n, table, status)
    end subroutine simple_set_profile_table

    !> Physical reference-coordinate metadata: lengths to metres, radial kind
    !> 1 = VMEC s, 2 = chart rho. Normalized TEST and SPECTRE are unsupported.
    subroutine simple_reference_info(length_to_m, radial_kind, status)
        use reference_coordinates, only: ref_coords
        use libneo_coordinates, only: chartmap_coordinate_system_t
        use simple_main, only: chartmap_cart_scale_to_m
        use magfie_sub, only: SPECTRE

        real(dp), intent(out) :: length_to_m
        integer, intent(out) :: radial_kind, status

        length_to_m = 0.01_dp
        radial_kind = 1
        status = 1
        if (isw_field_type == TEST .or. isw_field_type == SPECTRE) return
        if (.not. allocated(ref_coords)) return
        select type (cs => ref_coords)
        class is (chartmap_coordinate_system_t)
            length_to_m = chartmap_cart_scale_to_m
            if (length_to_m <= 0.0_dp) return
            radial_kind = 2
        class default
            continue
        end select
        status = 0
    end subroutine simple_reference_info

    !> Convert physical (R [m], phi [rad], Z [m]) to native reference coordinates.
    !> Return one inverse status per point; invalid results are NaN.
    subroutine simple_cylindrical_to_reference(n, cylindrical, reference, status)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
                                                                          ieee_quiet_nan
        use reference_coordinates, only: ref_coords

        integer, intent(in) :: n
        real(dp), intent(in) :: cylindrical(3, n)
        real(dp), intent(out) :: reference(3, n)
        integer, intent(out) :: status(n)
        real(dp) :: scale, point(3), result(3)
        integer :: i, radial_kind, info_status

        reference = ieee_value(0.0_dp, ieee_quiet_nan)
        status = -1
        call simple_reference_info(scale, radial_kind, info_status)
        if (info_status /= 0) return
        do i = 1, n
            if (.not. all(ieee_is_finite(cylindrical(:, i)))) cycle
            if (cylindrical(1, i) <= 0.0_dp) cycle
            point = cylindrical(:, i)
            point(1) = point(1)/scale
            point(3) = point(3)/scale
            call ref_coords%from_cyl(point, result, status(i))
            if (status(i) /= 0) cycle
            if (.not. all(ieee_is_finite(result))) then
                status(i) = -1
                cycle
            end if
            reference(:, i) = result
        end do
    end subroutine simple_cylindrical_to_reference

    !> Convert native reference coordinates to physical (R [m], phi [rad], Z [m]).
    subroutine simple_reference_to_cylindrical(n, reference, cylindrical, status)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
                                                                          ieee_quiet_nan
        use reference_coordinates, only: ref_coords

        integer, intent(in) :: n
        real(dp), intent(in) :: reference(3, n)
        real(dp), intent(out) :: cylindrical(3, n)
        integer, intent(out) :: status(n)
        real(dp) :: scale, point(3)
        integer :: i, radial_kind, info_status

        cylindrical = ieee_value(0.0_dp, ieee_quiet_nan)
        status = -1
        call simple_reference_info(scale, radial_kind, info_status)
        if (info_status /= 0) return
        do i = 1, n
            if (.not. all(ieee_is_finite(reference(:, i)))) cycle
            call ref_coords%evaluate_cyl(reference(:, i), point)
            if (.not. all(ieee_is_finite(point))) cycle
            point(1) = point(1)*scale
            point(3) = point(3)*scale
            cylindrical(:, i) = point
            status(i) = 0
        end do
    end subroutine simple_reference_to_cylindrical

    !> Restore all &config and &profiles values to their compiled defaults.
    subroutine simple_reset_config()
        call reset_config
        call reset_profiles
    end subroutine simple_reset_config

    !> Read the &config and &profiles groups of a simple.in file. Values not
    !> in the file keep their current setting. status is nonzero on error.
    subroutine simple_read_config(config_file, status, message)
        character(*), intent(in) :: config_file
        integer, intent(out) :: status
        character(len=256), intent(out) :: message

        integer :: unit

        message = ''
        open (newunit=unit, file=config_file, status='old', action='read', &
            iostat=status, iomsg=message)
        if (status /= 0) return
        read (unit, nml=config, iostat=status, iomsg=message)
        close (unit)
        if (status /= 0) return
        call read_profiles(config_file)
    end subroutine simple_read_config

    !> Set &config values from namelist text such as "&config ntestpart=8 /".
    !> status is nonzero for unknown names or malformed values.
    subroutine simple_set_config(namelist_text, status, message)
        character(*), intent(in) :: namelist_text
        integer, intent(out) :: status
        character(len=256), intent(out) :: message

        message = ''
        read (namelist_text, nml=config, iostat=status, iomsg=message)
    end subroutine simple_set_config

    !> Validate the configuration and load the field, as simple.x does before
    !> sampling.
    !> A Boozer chartmap input selects isw_field_type = 2 (Boozer) when neither
    !> isw_field_type nor integ_coords was set.
    subroutine simple_init()
        logical :: field_type_unset

        field_type_unset = integ_coords == -1000 .and. isw_field_type == CANFLUX
        call finalize_config
        if (field_type_unset .and. len_trim(field_input) > 0) then
            if (is_boozer_chartmap(field_input)) then
                isw_field_type = BOOZER
                integ_coords = BOOZER
            end if
        end if
        call validate_orbit_model_config
        call init_run(tracer)
        call init_sampling_field
        starting_line_stale = .false.
        active_field = FIELD_SAMPLING
        sampling_surface = sbeg(1)
    end subroutine simple_init

    !> Sample n start states as configured by startmode, sbeg and num_surf.
    subroutine simple_sample(n, z)
        integer, intent(in) :: n
        real(dp), intent(out) :: z(5, n)

        call allocate_particles(n)
        call use_sampling_field(needs_starting_line=.true.)
        call sample_start
        z = zstart(1:5, 1:n)
    end subroutine simple_sample

    !> Sample n start states on the flux surface s along the starting field line.
    subroutine simple_sample_surface(n, s, z)
        integer, intent(in) :: n
        real(dp), intent(in) :: s
        real(dp), intent(out) :: z(5, n)

        startmode = 1
        num_surf = 1
        sbeg(1) = s
        call simple_sample(n, z)
    end subroutine simple_sample_surface

    !> Sample n start states uniformly in the volume s_inner < s < s_outer.
    subroutine simple_sample_volume(n, s_inner, s_outer, z)
        integer, intent(in) :: n
        real(dp), intent(in) :: s_inner, s_outer
        real(dp), intent(out) :: z(5, n)

        startmode = 5
        num_surf = 2
        sbeg(1:2) = [s_inner, s_outer]
        call simple_sample(n, z)
    end subroutine simple_sample_volume

    !> Use the n states z as start states. When num_surf = 1 and all states
    !> share one radius, that radius becomes the starting surface, so bmin and
    !> bmax refer to it.
    subroutine simple_set_start(n, z)
        integer, intent(in) :: n
        real(dp), intent(in) :: z(5, n)

        call allocate_particles(n)
        zstart(1:5, 1:n) = z
        if (n > 0) then
            if (maxval(z(1, :)) <= minval(z(1, :))) call set_start_surface(z(1, 1))
        end if
    end subroutine simple_set_start

    !> Trace the start states to trace_time or loss, as simple.x does.
    subroutine simple_trace()
        call use_tracing_field
        call init_counters
        call trace_parallel(tracer)
    end subroutine simple_trace

    !> Trace the start states with orbit classification enabled.
    subroutine simple_classify()
        call use_tracing_field
        call init_counters
        call classify_parallel(tracer)
    end subroutine simple_classify

    !> Copy the per-particle results of the last trace or classification.
    !> times_lost is trace_time for confined, -1 for skipped and NaN for
    !> numerically failed particles.
    subroutine simple_get_results(n, z_end, loss_time, trap_parameter, &
            perp_invariant)
        integer, intent(in) :: n
        real(dp), intent(out) :: z_end(5, n), loss_time(n), trap_parameter(n), &
            perp_invariant(n)

        z_end = zend(1:5, 1:n)
        loss_time = times_lost(1:n)
        trap_parameter = trap_par(1:n)
        perp_invariant = perp_inv(1:n)
    end subroutine simple_get_results

    !> Copy the classification of the last simple_classify call. passing and
    !> lost are 0 or 1; classes holds the J_parallel, topological and fractal
    !> class (0 unclassified, 1 regular, 2 chaotic) per particle.
    subroutine simple_get_classification(n, passing, lost, classes)
        integer, intent(in) :: n
        integer, intent(out) :: passing(n), lost(n), classes(3, n)

        passing = merge(1, 0, class_passing(1:n))
        lost = merge(1, 0, class_lost(1:n))
        classes = iclass(1:3, 1:n)
    end subroutine simple_get_classification

    !> Write the simple.x output files for the last trace.
    subroutine simple_write_output()
        call write_output
    end subroutine simple_write_output

    !> Trace one orbit and return its ntimstep samples in reference coordinates.
    subroutine simple_trace_orbit(z0, nt, trajectory, times, z_end, loss_time)
        real(dp), intent(in) :: z0(5)
        integer, intent(in) :: nt
        real(dp), intent(out) :: trajectory(5, nt), times(nt), z_end(5)
        real(dp), intent(out) :: loss_time

        real(dp), allocatable :: integ_trajectory(:, :)
        integer :: it

        if (nt /= ntimstep) error stop 'simple_trace_orbit: nt must equal ntimstep'
        call simple_set_start(1, reshape(z0, [5, 1]))
        call use_tracing_field
        call init_counters
        allocate (integ_trajectory(5, nt))
        call trace_orbit(tracer, 1, integ_trajectory, times)
        do it = 1, nt
            call integ_to_ref(integ_trajectory(1:3, it), trajectory(1:3, it))
            trajectory(4:5, it) = integ_trajectory(4:5, it)
        end do
        z_end = zend(1:5, 1)
        loss_time = times_lost(1)
    end subroutine simple_trace_orbit

    !> Bounce/transit period, toroidal shift and frequencies of one orbit.
    !> values: period, period_std, delta_phi, delta_phi_std, omega_b, omega_phi.
    !> metadata: status, orbit_class, parallel_direction, n_periods, n_steps.
    subroutine simple_orbit_frequencies(z0, n_periods, max_steps, values, metadata)
        real(dp), intent(in) :: z0(5)
        integer, intent(in) :: n_periods, max_steps
        real(dp), intent(out) :: values(6)
        integer, intent(out) :: metadata(5)

        call set_start_surface(z0(1))
        call use_tracing_field
        call compute_canonical_frequencies_flat(tracer, z0, n_periods, max_steps, &
            values, metadata)
    end subroutine simple_orbit_frequencies

    !> Invariants of one state. values: h0, j_perp, p_phi, axis flux, outward
    !> sign, edge flux.
    subroutine simple_invariants_from_state(z0, values, status)
        real(dp), intent(in) :: z0(5)
        real(dp), intent(out) :: values(6)
        integer, intent(out) :: status

        call use_tracing_field
        call invariants_from_state_flat(tracer, z0, values, status)
    end subroutine simple_invariants_from_state

    !> Start states on the axisymmetric cut matching (h0, j_perp, p_phi),
    !> ordered by major radius. metadata: sigma, section branch, section kind.
    subroutine simple_states_from_invariants(invariants, max_solutions, states, &
            metadata, residuals, cylindrical, n_solutions, status)
        real(dp), intent(in) :: invariants(3)
        integer, intent(in) :: max_solutions
        real(dp), intent(out) :: states(5, max_solutions)
        integer, intent(out) :: metadata(3, max_solutions)
        real(dp), intent(out) :: residuals(max_solutions)
        real(dp), intent(out) :: cylindrical(3, max_solutions)
        integer, intent(out) :: n_solutions, status

        call use_tracing_field
        call states_from_invariants_flat(tracer, invariants, max_solutions, &
            states, metadata, residuals, cylindrical, n_solutions, status)
    end subroutine simple_states_from_invariants

    !> Convert POTATO invariants (h0, j_perp, psi_star) to SIMPLE's convention.
    subroutine simple_potato_invariants(potato, psi_axis, psi_edge, v0_ratio, &
            invariants)
        real(dp), intent(in) :: potato(3), psi_axis, psi_edge, v0_ratio
        real(dp), intent(out) :: invariants(3)

        call potato_invariants_to_simple_flat(potato, psi_axis, psi_edge, &
            v0_ratio, invariants)
    end subroutine simple_potato_invariants

    !> Trace one orbit to its next tip (cut 0), toroidal-period cut (1) or
    !> either (-1). cut_state holds the state and the parallel invariant.
    subroutine simple_trace_to_cut(z0, cut, max_events, cut_state, cut_type, &
            status)
        real(dp), intent(in) :: z0(5)
        integer, intent(in) :: cut, max_events
        real(dp), intent(out) :: cut_state(6)
        integer, intent(out) :: cut_type, status

        call set_start_surface(z0(1))
        call use_tracing_field
        call trace_to_cut_flat(tracer, z0, cut, max_events, cut_state, cut_type, &
            status)
    end subroutine simple_trace_to_cut

    subroutine allocate_particles(n)
        integer, intent(in) :: n

        ntestpart = n
        call reallocate_arrays
        starting_line_stale = .true.
    end subroutine allocate_particles

    subroutine set_start_surface(s)
        real(dp), intent(in) :: s

        if (isw_field_type == TEST) then
            ! Analytic tokamak: B = 1 - r cos(theta), so B ranges over 1 -+ r.
            bmod00 = 1.0_dp
            bmin = 1.0_dp - s
            bmax = 1.0_dp + s
        else if (num_surf == 1) then
            sbeg(1) = s
            call use_sampling_field
        end if
    end subroutine set_start_surface

    subroutine use_sampling_field(needs_starting_line)
        ! The starting field line depends on sbeg(1); integrate it only when the
        ! sampling field is not active on that surface, or when the caller
        ! samples and the line arrays were reset by reallocate_arrays.
        logical, intent(in), optional :: needs_starting_line
        logical :: refill

        if (active_field == FIELD_NONE) error stop 'simple_api: call simple_init first'
        refill = .false.
        if (present(needs_starting_line)) refill = needs_starting_line .and. &
                                                   starting_line_stale
        if (active_field == FIELD_SAMPLING .and. .not. refill) then
            if (abs(sbeg(1) - sampling_surface) <= 0.0_dp) return
        end if
        call init_sampling_field
        starting_line_stale = .false.
        active_field = FIELD_SAMPLING
        sampling_surface = sbeg(1)
    end subroutine use_sampling_field

    subroutine use_tracing_field
        if (active_field == FIELD_NONE) error stop 'simple_api: call simple_init first'
        if (active_field /= FIELD_TRACING) call init_tracing_field
        active_field = FIELD_TRACING
    end subroutine use_tracing_field

end module simple_api
