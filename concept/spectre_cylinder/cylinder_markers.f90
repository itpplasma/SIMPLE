program cylinder_markers
    ! Production GC/FO hybrid (integmode = 3, crossing_level = 3) for markers
    ! launched inside the FO entry band of interface 1. Each marker records its
    ! entry particle, whether the particle reaches the other volume within
    ! NGYRO gyro-periods, and mu and the guiding centre at the first microstep
    ! end on the far side. Usage: cylinder_markers.x <h5> <nmarker> <out>
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple_api, only: simple_read_config, simple_init
    use simple, only: init_sympl
    use params, only: relerr, integmode, crossing_level
    use magfie_sub, only: init_magfie, SPECTRE, spectre_field, M_TO_CM, &
                          TESLA_TO_GAUSS
    use field_can_mod, only: field_can_t, ref_to_integ
    use field_can_spectre, only: spectre_mvol, set_spectre_volume_lock
    use orbit_symplectic_base, only: symplectic_integrator_t
    use spectre_fo_field, only: spectre_fo_field_local, spectre_fo_local_t
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_enter, &
                                 spectre_fo_radial_larmor, SPECTRE_FO_OK
    use spectre_sympl_orbit, only: sympl_spectre_state_t, sympl_spectre_reset, &
                                   orbit_microstep_sympl_spectre, SYMPL_SPECTRE_OK
    use parmot_mod, only: ro0
    use util, only: sqrt2, twopi
    implicit none

    integer, parameter :: NPER = 32, NGYRO = 3
    real(dp), parameter :: BAND = 1.4_dp, LAMBDA_MAX = 0.9_dp
    character(len=1024) :: fixture, outfile, arg
    character(len=256) :: message
    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f
    type(sympl_spectre_state_t) :: state
    type(spectre_fo_state_t) :: entry
    type(spectre_fo_local_t) :: loc
    real(dp) :: y(5), z(5), rnd(4), rho_r, dt, r0, b_ref, fraction
    real(dp) :: r_gc, r_if, b_gc, mu0, mu1, r_x, r_if1, vpar, xgc(3), d0
    integer :: nmark, i, k, status, unit, owner, side0, crossed, kcross

    call get_command_argument(1, fixture)
    call get_command_argument(2, arg)
    read (arg, *) nmark
    call get_command_argument(3, outfile)
    open (newunit=unit, file='simple.in', status='replace')
    write (unit, '(a)') '&config'
    write (unit, '(a)') "field_input='"//trim(fixture)//"'"
    write (unit, '(a)') 'integ_coords=6, integmode=3, crossing_level=3'
    write (unit, '(a)') 'ntestpart=1, ntimstep=65, npoiper2=1024'
    write (unit, '(a)') 'facE_al=1, trace_time=1d-5, relerr=1d-12 /'
    close (unit)
    call simple_read_config('simple.in', status, message)
    if (status /= 0) error stop 'configuration'
    call simple_init
    call init_magfie(SPECTRE)

    r0 = spectre_field%data%Rbc(findloc(spectre_field%data%im, 0, dim=1), &
                                spectre_mvol)
    call spectre_fo_field_local([1.0_dp - 1d-9, 0.0_dp, 0.0_dp], 1, loc)
    b_ref = loc%bmod
    ! Microstep: 1/NPER of the gyro-period at the sheet, in path length (cm).
    dt = twopi*ro0/b_ref/NPER

    call random_seed(put=[(12345 + 7*i, i=1, 64)])
    open (newunit=unit, file=trim(outfile), status='replace', action='write')
    write (unit, '(A,ES24.16,A,ES24.16,A,ES24.16)') '# ro0_cmG ', ro0, &
        ' dt_cm ', dt, ' R0_m ', r0
    write (unit, '(A)') '# i side_gc lambda theta zeta d_rhor r_gc_m '// &
        'd_gc_m B_gc_T mu0 x0_cm(3) v0(3) side0 crossed t_cross_cm mu1 '// &
        'r_after_m d_after_m status d0_m'
    do i = 1, nmark
        call random_number(rnd)
        y = [1.0_dp, twopi*rnd(1), twopi*rnd(2), 1.0_dp, &
             LAMBDA_MAX*(2*rnd(3) - 1)]
        owner = merge(1, 2, rnd(4) < 0.5_dp)
        call spectre_fo_radial_larmor(y, owner, ro0/sqrt2, rho_r, status)
        if (status /= SPECTRE_FO_OK) error stop 'radial gyroradius'
        ! Uniform in |d| < BAND rho_r on the owner's side of the sheet.
        call random_number(rnd(1))
        y(1) = 1.0_dp + merge(-1.0_dp, 1.0_dp, owner == 1)*BAND*rnd(1)*rho_r
        call spectre_fo_field_local(y(1:3), owner, loc)
        r_gc = radius(loc%x)
        b_gc = loc%bmod/TESLA_TO_GAUSS
        mu0 = (1 - y(5)**2)/b_gc
        call spectre_fo_field_local([1.0_dp, y(2), y(3)], owner, loc)
        r_if = radius(loc%x)
        ! The production entry draws the same hashed gyrophase.
        call spectre_fo_enter(entry, y, owner, 1, ro0/sqrt2, status)
        if (status /= SPECTRE_FO_OK) error stop 'entry'
        side0 = entry%owner
        call spectre_fo_field_local([1.0_dp, theta_of(entry%x), entry%u(3)], &
                                    side0, loc)
        d0 = radius(entry%x) - radius(loc%x)

        call set_spectre_volume_lock(0)
        call ref_to_integ(y(1:3), z(1:3))
        z(4:5) = y(4:5)
        call init_sympl(si, f, z, dt, dt, relerr, integmode)
        call sympl_spectre_reset(state, si, spectre_mvol, integmode, &
                                 crossing_level)
        crossed = 0
        kcross = 0
        mu1 = 0; r_x = 0; r_if1 = 0
        do k = 1, NGYRO*NPER
            call orbit_microstep_sympl_spectre(state, si, f, status, fraction)
            if (status /= SYMPL_SPECTRE_OK) exit
            if (.not. state%fo%active) then
                status = -1
                exit
            end if
            ! The hybrid entered at k = 1 from the same state as entry.
            if (k == 1 .and. (abs(state%fo%entry_rho - entry%entry_rho) > &
                            1d-12 .or. abs(state%fo%entry_pzeta/entry%entry_pzeta - 1) &
                              > 1d-9)) error stop 'entry mismatch'
            if (state%fo%owner /= side0) then
                crossed = 1
                kcross = k
                call spectre_fo_field_local(state%fo%u, state%fo%owner, loc)
                vpar = dot_product(state%fo%v, loc%b)
                mu1 = (dot_product(state%fo%v, state%fo%v) - vpar**2)/2/ &
                      (loc%bmod/TESLA_TO_GAUSS)
                xgc = state%fo%x - ro0/sqrt2/loc%bmod*cross(loc%b, state%fo%v)
                r_x = radius(xgc)
                call spectre_fo_field_local([1.0_dp, theta_of(xgc), &
                                             state%fo%u(3)], state%fo%owner, loc)
                r_if1 = radius(loc%x)
                exit
            end if
        end do
        write (unit, '(I6,I3,8ES16.8,6ES18.10,2I3,4ES16.8,I4,ES16.8)') &
            i, owner, y(5), y(2), y(3), (y(1) - 1)/rho_r, r_gc, r_gc - r_if, &
            b_gc, mu0, entry%x, entry%v, side0, crossed, kcross*dt, mu1, r_x, &
            r_x - r_if1, status, d0
    end do
    close (unit)

contains

    real(dp) function radius(x) result(r)
        real(dp), intent(in) :: x(3)
        r = hypot(hypot(x(1), x(2))/M_TO_CM - r0, x(3)/M_TO_CM)
    end function radius

    real(dp) function theta_of(x) result(t)
        ! Right-handed angle about e_phi: theta = 0 outboard, Z < 0 at pi/2.
        real(dp), intent(in) :: x(3)
        t = modulo(atan2(-x(3)/M_TO_CM, hypot(x(1), x(2))/M_TO_CM - r0), twopi)
    end function theta_of

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function cross
end program cylinder_markers
