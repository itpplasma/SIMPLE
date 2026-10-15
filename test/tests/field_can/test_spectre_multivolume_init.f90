program test_spectre_multivolume_init
    ! Production initialization must permit GC motion in the outer volume,
    ! without hiding a Newton-domain failure behind bulk full-orbit recovery.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use simple_api, only: simple_read_config, simple_init
    use simple, only: init_sympl
    use params, only: dtaumin, relerr, integmode, crossing_level
    use field_can_mod, only: field_can_t, ref_to_integ
    use field_can_spectre, only: spectre_mvol
    use magfie_sub, only: init_magfie, SPECTRE
    use orbit_symplectic_base, only: symplectic_integrator_t, sympl_rmax
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_enter, &
        spectre_fo_advance_until_exit, SPECTRE_FO_OK, SPECTRE_FO_LOSS
    use parmot_mod, only: ro0
    use spectre_sympl_orbit, only: sympl_spectre_state_t, sympl_spectre_reset, &
        orbit_microstep_sympl_spectre, sympl_landing_stats_reset, &
        sympl_landing_stats, sympl_fo_stats, SYMPL_SPECTRE_OK, SYMPL_SPECTRE_LOSS
    use util, only: twopi, sqrt2

    implicit none

    integer, parameter :: NSTEP = 10000
    character(len=1024) :: h5file
    character(len=256) :: message
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f
    type(sympl_spectre_state_t) :: state
    type(spectre_fo_state_t) :: reference
    real(dp) :: z(5), x(3), initial_angles(2), advance, t_frac, max_resid
    real(dp) :: y(5), used
    integer :: unit, status, k, completed, entries, exits, losses, failures
    integer :: counts(5), landings, stops
    integer :: reference_status, reference_owner
    logical :: exited

    call get_command_argument(1, h5file)
    if (len_trim(h5file) == 0) error stop 'missing tok2vol.h5 fixture'
    open (newunit=unit, file='simple.in', status='replace', action='write')
    write (unit, '(A)') '&config'
    write (unit, '(A)') " field_input = '"//trim(h5file)//"'"
    write (unit, '(A)') ' integ_coords = 6, integmode = 3, crossing_level = 3'
    write (unit, '(A)') ' ntestpart = 1, ntimstep = 100, npoiper2 = 1024'
    write (unit, '(A)') ' facE_al = 1000d0, trace_time = 1d-3, relerr = 1d-13'
    write (unit, '(A)') ' spectre_ncon_phi = 32'
    write (unit, '(A)') '/'
    close (unit)
    call simple_read_config('simple.in', status, message)
    if (status /= 0) then
        print *, trim(message)
        error stop 'configuration failed'
    end if
    call simple_init
    call init_magfie(SPECTRE)
    if (spectre_mvol /= 2) error stop 'expected a two-volume fixture'
    print '(A,F5.1)', 'production sympl_rmax = ', sympl_rmax

    x = [1.5_dp, 0.3_dp, 0.0_dp]
    call ref_to_integ(x, z(1:3))
    z(4:5) = [1.0_dp, 0.8_dp]
    call init_sympl(si, f, z, dtaumin, dtaumin, relerr, integmode)
    call sympl_spectre_reset(state, si, spectre_mvol, integmode, crossing_level)
    call sympl_landing_stats_reset
    initial_angles = si%z(2:3)
    completed = 0
    do k = 1, NSTEP
        call orbit_microstep_sympl_spectre(state, si, f, status, t_frac)
        if (status /= SYMPL_SPECTRE_OK) exit
        if (.not. all(ieee_is_finite(si%z))) error stop 'nonfinite GC state'
        if (abs(t_frac - 1.0_dp) > 1d-12) error stop 'incomplete microstep'
        if (si%z(1) <= 1.1_dp .or. si%z(1) >= 1.9_dp) &
            error stop 'marker left the bulk of volume 2'
        completed = completed + 1
        call sympl_fo_stats(entries, exits, losses, failures, counts)
        if (entries /= 0 .or. state%warning_hold_count /= 0) exit
    end do
    call sympl_fo_stats(entries, exits, losses, failures, counts)
    call sympl_landing_stats(landings, max_resid, stops)
    advance = maxval(abs(si%z(2:3) - initial_angles))
    print '(A,I0,A,I0,A,I0,A,I0)', 'GC steps = ', completed, &
        ', FO entries = ', entries, ', stops = ', stops, &
        ', warning holds = ', state%warning_hold_count
    print '(A,ES12.4)', 'angular advance = ', advance
    if (completed /= NSTEP .or. entries /= 0 .or. stops /= 0 .or. &
        state%warning_hold_count /= 0 .or. advance < 1.0_dp) &
        error stop 'outer-volume GC did not advance normally'
    print *, 'production multi-volume GC PASS'

    ! The particle, rather than its GC, must hit the absorbing edge. Compare
    ! the hybrid against a directly launched FO (whose edge event is checked
    ! independently by test_spectre_fo_events). The centre remains inside.
    x = [2.0_dp - 1d-7, 3.0_dp*twopi/4.0_dp, 0.0_dp]
    call ref_to_integ(x, z(1:3))
    z(4:5) = [1.0_dp, 0.8_dp]
    call init_sympl(si, f, z, dtaumin, dtaumin, relerr, integmode)
    call sympl_spectre_reset(state, si, spectre_mvol, integmode, crossing_level)
    call sympl_landing_stats_reset
    y = [x, z(4:5)]
    call spectre_fo_enter(reference, y, spectre_mvol, spectre_mvol, &
        ro0/sqrt2, reference_status)
    if (reference_status /= SPECTRE_FO_OK) error stop 'reference edge launch'
    reference%exit_rho = 1.0e6_dp
    do k = 1, 100
        call spectre_fo_advance_until_exit(reference, state%dt_std, ro0/sqrt2, &
            used, y, reference_owner, exited, reference_status)
        call orbit_microstep_sympl_spectre(state, si, f, status, t_frac)
        if (.not. state%fo%active) error stop 'outer edge did not enter FO'
        if (maxval(abs(state%fo%x - reference%x)) > 1d-10 .or. &
            maxval(abs(state%fo%v - reference%v)) > 1d-12) &
            error stop 'hybrid wall orbit differs from direct FO'
        if (abs(t_frac - used/state%dt_std) > 1d-12) &
            error stop 'hybrid wall event time differs from direct FO'
        if (status /= SYMPL_SPECTRE_OK) exit
    end do
    print '(A,I0,A,F14.10,A,ES12.4)', 'outer-edge status = ', status, &
        ', radius = ', si%z(1), ', time fraction = ', t_frac
    if (status /= SYMPL_SPECTRE_LOSS) error stop 'outer marker not lost'
    if (reference_status /= SPECTRE_FO_LOSS) error stop 'reference not lost'
    if (abs(state%fo%u(1) - real(spectre_mvol, dp)) > 1d-8) &
        error stop 'particle loss not localized at Mvol'
    if (si%z(1) >= real(spectre_mvol, dp)) &
        error stop 'missing finite-Larmor-radius wall-loss witness'
    if (.not. ieee_is_finite(t_frac)) error stop 'invalid loss fraction'
    if (t_frac <= 0.0_dp .or. t_frac > 1.0_dp) error stop 'invalid loss time'
    call sympl_fo_stats(entries, exits, losses, failures, counts)
    call sympl_landing_stats(landings, max_resid, stops)
    if (entries /= 1 .or. losses /= 1 .or. failures /= 0 .or. &
        stops /= 0 .or. landings /= 0 .or. state%warning_hold_count /= 0) &
        error stop 'outer loss used numerical recovery'
    print *, 'physical outer-edge loss PASS'
end program test_spectre_multivolume_init
