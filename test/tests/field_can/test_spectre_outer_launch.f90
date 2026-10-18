program test_spectre_outer_launch
    ! A GC inside the absorbing wall can represent a particle already outside.
    ! Check geometry independently with the raw boundary and its outward normal.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple_api, only: simple_read_config, simple_init
    use simple, only: init_sympl
    use params, only: dtaumin, relerr, integmode, crossing_level
    use magfie_sub, only: init_magfie, SPECTRE, spectre_field, M_TO_CM
    use field_can_mod, only: field_can_t, ref_to_integ
    use field_can_spectre, only: spectre_mvol, set_spectre_volume_lock
    use orbit_symplectic_base, only: symplectic_integrator_t
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_enter, &
        spectre_fo_radial_larmor, SPECTRE_FO_OK, SPECTRE_FO_LOSS
    use spectre_sympl_orbit, only: sympl_spectre_state_t, sympl_spectre_reset, &
        orbit_microstep_sympl_spectre, sympl_landing_stats_reset, sympl_fo_stats, &
        SYMPL_SPECTRE_LOSS
    use parmot_mod, only: ro0
    use util, only: sqrt2, twopi
    implicit none

    character(len=1024) :: fixture
    character(len=256) :: message
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f
    type(sympl_spectre_state_t) :: state
    type(spectre_fo_state_t) :: reference
    real(dp) :: y(5), z(5), rho_l, fraction, wall_u(3), wall_x(3), e(3, 3), n(3)
    integer :: unit, status, j, outside, inside, entries, exits, losses, failures
    integer :: counts(5)

    call get_command_argument(1, fixture)
    open (newunit=unit, file='simple.in', status='replace')
    write (unit, '(a)') '&config'
    write (unit, '(a)') "field_input='"//trim(fixture)//"'"
    write (unit, '(a)') 'integ_coords=6, integmode=3, crossing_level=3'
    write (unit, '(a)') 'ntestpart=1, ntimstep=65, npoiper2=1024'
    write (unit, '(a)') 'facE_al=1000, trace_time=2d-3, relerr=1d-12 /'
    close (unit)
    call simple_read_config('simple.in', status, message)
    if (status /= 0) error stop 'outer-launch configuration'
    call simple_init
    call init_magfie(SPECTRE)
    outside = 0
    inside = 0
    do j = 1, 16
        call set_spectre_volume_lock(0)
        y = [real(spectre_mvol, dp) - 1d-7, twopi*(j - 1)/16, &
            0.37_dp, 1.0_dp, 0.4_dp]
        call spectre_fo_radial_larmor(y, spectre_mvol, ro0/sqrt2, rho_l, status)
        if (status /= SPECTRE_FO_OK) error stop 'outer-launch gyroradius'
        y(1) = real(spectre_mvol, dp) - 0.25_dp*rho_l
        call spectre_fo_enter(reference, y, spectre_mvol, spectre_mvol, &
            ro0/sqrt2, status)
        if (status == SPECTRE_FO_OK) then
            inside = inside + 1
            cycle
        end if
        if (status /= SPECTRE_FO_LOSS) error stop 'outside launch misclassified'
        wall_u = [real(spectre_mvol, dp), reference%u(2:3)]
        call spectre_field%coords%evaluate_cart(wall_u, wall_x)
        call spectre_field%coords%covariant_basis(wall_u, e)
        n = cross_product(e(:, 2), e(:, 3))
        n = n*sign(1.0_dp, dot_product(n, e(:, 1)))/norm2(n)
        if (dot_product(reference%x - wall_x*M_TO_CM, n) <= 1d-5) &
            error stop 'raw boundary does not confirm exterior particle'
        outside = outside + 1
        call ref_to_integ(y(1:3), z(1:3))
        z(4:5) = y(4:5)
        call init_sympl(si, f, z, dtaumin, dtaumin, relerr, integmode)
        call sympl_spectre_reset(state, si, spectre_mvol, integmode, crossing_level)
        call sympl_landing_stats_reset
        call orbit_microstep_sympl_spectre(state, si, f, status, fraction)
        if (status /= SYMPL_SPECTRE_LOSS .or. abs(fraction) > 1d-12) &
            error stop 'outside particle is not lost immediately'
        if (maxval(abs(state%fo%x - reference%x)) > 1d-7) &
            error stop 'wrong particle at initial wall loss'
        call sympl_fo_stats(entries, exits, losses, failures, counts)
        if (entries /= 1 .or. losses /= 1 .or. failures /= 0 .or. &
            state%warning_hold_count /= 0) &
            error stop 'physical loss counted as failure'
    end do
    if (outside == 0 .or. inside == 0) error stop 'missing wall-launch geometry'
    print '(a,2i4)', 'inside/outside particle launches: ', inside, outside

contains

    pure function cross_product(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function cross_product
end program test_spectre_outer_launch
