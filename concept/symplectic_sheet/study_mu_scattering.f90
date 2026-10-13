program study_mu_scattering
    !> Statistics of the magnetic-moment change of markers that cross the sheet
    !> thousands of times, for one scheme per run.
    !>
    !> usage: study_mu_scattering.x SCHEME H WIDTH A_OUT P_PHI [G_MIN L]
    !>   SCHEME  rk45 | midpoint_region | naive | landing | poincare
    !>   H       step (rk45: tolerance)
    !>   WIDTH   smooth layer width (0: sharp sheet; poincare needs > 0)
    !>   A_OUT   outer psi coefficient (inner is 1): sets the direction jump
    !>   G_MIN L optional Poincare time-scale floor and length (default 0.1,
    !>           2 WIDTH)
    !>   P_PHI   canonical toroidal momentum; v_phi at the sheet is
    !>           (P_PHI - psi(r0))/R and sets the force jump v_phi x [[B_p]]
    !>
    !> Each marker starts next to the sheet; mu is time-averaged over the first
    !> and last T_AVG of T_END. Printed are the rms and mean of the relative
    !> change over the ensemble, the crossings per marker and the steps per
    !> unit time.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sheet_field, only: field_params_t, region_of, bfield, psi_value, &
        hamiltonian
    use fo_sheet, only: sheet_step, step_stats_t, SCHEME_NAIVE, &
        SCHEME_MIDPOINT_REGION, SCHEME_LANDING, poincare_step, poincare_params_t
    use rk_reference, only: rk_advance
    use sheet_setup, only: initial_state
    implicit none

    integer, parameter :: NMARKER = 256
    real(dp), parameter :: T_END = 4.0e4_dp, T_AVG = 200.0_dp, SPEED = 0.05_dp
    real(dp), parameter :: RK_SAMPLE = 0.1_dp

    type(field_params_t) :: par
    type(poincare_params_t) :: pp
    character(32) :: scheme, arg
    real(dp) :: h, mu0(NMARKER), mu1(NMARKER), d(NMARKER)
    integer :: im, crossings(NMARKER)
    integer(8) :: steps(NMARKER)

    call get_command_argument(1, scheme)
    call get_command_argument(2, arg)
    read (arg, *) h
    call get_command_argument(3, arg)
    read (arg, *) par%width
    call get_command_argument(4, arg)
    read (arg, *) par%a_out
    call get_command_argument(5, arg)
    read (arg, *) par%p_phi
    pp%g_min = 0.1_dp
    pp%length = 2.0_dp*par%width
    if (command_argument_count() >= 7) then
        call get_command_argument(6, arg)
        read (arg, *) pp%g_min
        call get_command_argument(7, arg)
        read (arg, *) pp%length
    end if

    !$omp parallel do default(shared) private(im) schedule(dynamic)
    do im = 1, NMARKER
        call run_marker(im, mu0(im), mu1(im), crossings(im), steps(im))
    end do
    !$omp end parallel do

    d = mu1/mu0 - 1.0_dp
    print '(A16,A,F7.4,A,ES8.1,A,F4.1,A,F5.2,A,ES10.3,A,ES10.3,A,I6,A,F8.1)', &
        trim(scheme), ' h ', h, ' w ', par%width, ' a_out ', par%a_out, &
        ' p_phi ', par%p_phi, '  rms ', sqrt(sum(d**2)/NMARKER), '  mean ', sum(d)/NMARKER, &
        '  crossings ', sum(crossings)/NMARKER, '  steps/time ', &
        real(sum(steps), dp)/(NMARKER*T_END)

contains

    subroutine run_marker(im, m0, m1, ncross, nstep)
        integer, intent(in) :: im
        real(dp), intent(out) :: m0, m1
        integer, intent(out) :: ncross
        integer(8), intent(out) :: nstep

        real(dp) :: z(4), z_new(4), t, dt, h_try, e0
        integer :: region
        type(step_stats_t) :: stats
        logical :: ok

        call initial_state(par, SPEED, mod(im - 1, 8) + 1, z)
        z(3) = z(3)*(1.0_dp + 1.0e-3_dp*real(im, dp)/NMARKER)
        region = region_of(par, z(1:2))
        e0 = hamiltonian(par, region, z)
        m0 = 0.0_dp
        m1 = 0.0_dp
        t = 0.0_dp
        ncross = 0
        nstep = 0
        h_try = 0.01_dp
        do while (t < T_END)
            select case (trim(scheme))
            case ('rk45')
                dt = RK_SAMPLE
                call rk_advance(par, region, z, dt, h, h_try, ncross)
            case ('midpoint_region')
                dt = h
                call sheet_step(par, SCHEME_MIDPOINT_REGION, region, z, h, &
                    z_new, stats)
                z = z_new
            case ('naive')
                dt = h
                call sheet_step(par, SCHEME_NAIVE, region, z, h, z_new, stats)
                z = z_new
            case ('landing')
                dt = h
                call sheet_step(par, SCHEME_LANDING, region, z, h, z_new, stats)
                z = z_new
            case ('poincare')
                call poincare_step(par, pp, z, e0, h, z_new, dt, ok)
                if (.not. ok) stats%failures = stats%failures + 1
                if (region_of(par, z_new(1:2)) /= region) &
                    stats%crossings = stats%crossings + 1
                region = region_of(par, z_new(1:2))
                z = z_new
            case default
                error stop 'unknown scheme'
            end select
            t = t + dt
            nstep = nstep + 1
            if (t <= T_AVG) m0 = m0 + mu_of(z, region)*dt
            if (t > T_END - T_AVG) m1 = m1 + mu_of(z, region)*dt
        end do
        if (trim(scheme) /= 'rk45') ncross = stats%crossings
    end subroutine run_marker

    real(dp) function mu_of(z, region) result(mu)
        real(dp), intent(in) :: z(4)
        integer, intent(in) :: region

        real(dp) :: b(3), v(3), bm, vpar

        call bfield(par, z(1:2), b)
        v = [z(3), (par%p_phi - psi_value(par, region, z(1:2)))/z(1), &
             z(4) + par%f_tor*log(z(1))]
        bm = norm2(b)
        vpar = dot_product(v, b)/bm
        mu = (dot_product(v, v) - vpar**2)/(2.0_dp*bm)
    end function mu_of

end program study_mu_scattering
