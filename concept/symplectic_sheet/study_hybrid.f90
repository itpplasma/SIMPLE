program study_hybrid
    !> GC in the volumes, FO near the sheet: statistics against FO everywhere.
    !>
    !> usage: study_hybrid.x FO LIFT SPEED D_IN D_OUT T_END A_IN A_OUT [H]
    !>   FO     rk45 (sharp sheet, tol 1e-10) | gauss2 (Gauss-Legendre 2 with
    !>          Poincare refinement on an analytic layer of width SPEED/20)
    !>   LIFT   gyrophase | pphi (see gc_hybrid)
    !>   SPEED  particle speed (gyroradius SPEED at |B| = 1)
    !>   D_IN   switch GC -> FO when the GC is closer to the sheet (units of the
    !>          gyroradius)
    !>   D_OUT  switch FO -> GC when the GC of the FO state is farther; a huge
    !>          value keeps FO everywhere (the reference)
    !>   T_END  physical time
    !>   A_IN, A_OUT  psi coefficients inside and outside the sheet: the
    !>          poloidal field (safety factor) and the direction jump
    !>   H      fixed step of gauss2 (default 0.1)
    !>
    !> Markers start in FO next to the sheet with p_phi such that
    !> v_phi = 0.4 SPEED at the sheet.
    !> Printed: rms and mean relative change of the time-averaged mu, rms
    !> change of the FO p_phi (exactly conserved physically), rms relative
    !> energy error, switches and FO fraction of the time, wall time.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sheet_field, only: field_params_t, region_of, hamiltonian
    use fo_sheet, only: poincare_gauss2_step, poincare_params_t
    use rk_reference, only: rk_advance
    use sheet_setup, only: initial_state
    use gc_hybrid, only: gc_state_t, fo_to_gc, gc_to_fo, gc_advance, &
        LIFT_GYROPHASE, LIFT_PPHI
    implicit none

    integer, parameter :: NMARKER = 256
    real(dp), parameter :: T_AVG = 200.0_dp, RK_TOL = 1.0e-10_dp
    real(dp), parameter :: GC_TOL = 1.0e-10_dp, FO_SAMPLE = 0.1_dp

    type(field_params_t) :: par
    type(poincare_params_t) :: pp
    character(32) :: fo_scheme, lift_name, arg
    real(dp) :: speed, d_in, d_out, t_end, h
    real(dp) :: mu0(NMARKER), mu1(NMARKER), dpphi(NMARKER), de(NMARKER)
    real(dp) :: fo_frac(NMARKER), dmu(NMARKER)
    integer :: switches(NMARKER), lift_fail(NMARKER), lift, im
    integer(8) :: clock0, clock1, rate

    call get_command_argument(1, fo_scheme)
    call get_command_argument(2, lift_name)
    call get_command_argument(3, arg)
    read (arg, *) speed
    call get_command_argument(4, arg)
    read (arg, *) d_in
    call get_command_argument(5, arg)
    read (arg, *) d_out
    call get_command_argument(6, arg)
    read (arg, *) t_end
    call get_command_argument(7, arg)
    read (arg, *) par%a_in
    call get_command_argument(8, arg)
    read (arg, *) par%a_out
    h = 0.1_dp
    if (command_argument_count() >= 9) then
        call get_command_argument(9, arg)
        read (arg, *) h
    end if
    lift = merge(LIFT_PPHI, LIFT_GYROPHASE, trim(lift_name) == 'pphi')
    d_in = d_in*speed
    d_out = d_out*speed
    par%p_phi = par%a_in*par%r_sheet**2 + 0.4_dp*speed*par%r_major
    if (trim(fo_scheme) == 'gauss2') then
        par%width = speed/20.0_dp
        par%analytic = .true.
        pp%g_min = 0.1_dp
        pp%length = 0.01_dp*speed/0.05_dp
    end if

    call system_clock(clock0, rate)
    !$omp parallel do default(shared) private(im) schedule(dynamic)
    do im = 1, NMARKER
        call run_marker(im, mu0(im), mu1(im), dpphi(im), de(im), fo_frac(im), &
            switches(im), lift_fail(im))
    end do
    !$omp end parallel do
    call system_clock(clock1)

    dmu = mu1/mu0 - 1.0_dp
    print '(A8,A10,A,F5.3,A,2ES9.2,A,ES10.3,A,ES10.3,A,ES9.2,A,ES9.2,A,F9.1,A,F5.3,A,I5,A,F8.1)', &
        trim(fo_scheme), trim(lift_name), ' v ', speed, ' d_in,out ', d_in, &
        d_out, '  rms dmu ', sqrt(sum(dmu**2)/NMARKER), '  mean ', &
        sum(dmu)/NMARKER, '  rms dpphi ', sqrt(sum(dpphi**2)/NMARKER), &
        '  rms dH/H ', sqrt(sum(de**2)/NMARKER), '  switches ', &
        real(sum(switches), dp)/NMARKER, '  FO frac ', sum(fo_frac)/NMARKER, &
        '  lift fails ', sum(lift_fail), '  wall s ', &
        real(clock1 - clock0, dp)/real(rate, dp)

contains

    subroutine run_marker(im, m0, m1, dp_phi, d_e, frac, nswitch, nfail)
        integer, intent(in) :: im
        real(dp), intent(out) :: m0, m1, dp_phi, d_e, frac
        integer, intent(out) :: nswitch, nfail

        type(field_params_t) :: pm
        type(gc_state_t) :: gc
        real(dp) :: z(4), z_new(4), t, dt, e0, pphi0, h_try, t_fo, dt_gc, pphi
        real(dp) :: mu_now
        integer :: region, ncross
        logical :: in_fo, entered, ok

        pm = par
        call initial_state(pm, speed, mod(im - 1, 8) + 1, z)
        z(3) = z(3)*(1.0_dp + 1.0e-3_dp*real(im, dp)/NMARKER)
        region = region_of(pm, z(1:2))
        e0 = hamiltonian(pm, region, z)
        pphi0 = pm%p_phi
        in_fo = .true.
        t = 0.0_dp
        t_fo = 0.0_dp
        m0 = 0.0_dp
        m1 = 0.0_dp
        nswitch = 0
        nfail = 0
        ncross = 0
        h_try = 0.01_dp
        do while (t < t_end)
            if (in_fo) then
                if (trim(fo_scheme) == 'rk45') then
                    dt = FO_SAMPLE
                    call rk_advance(pm, region, z, dt, RK_TOL, h_try, ncross)
                else
                    call poincare_gauss2_step(pm, pp, z, e0, h, z_new, dt, ok)
                    z = z_new
                    region = region_of(pm, z(1:2))
                end if
                t_fo = t_fo + dt
                gc = fo_to_gc(pm, z)
                mu_now = gc%mu
                if (abs(sheet_s(gc)) > d_out) then
                    in_fo = .false.
                    nswitch = nswitch + 1
                end if
            else
                call gc_advance(pm, gc, min(FO_SAMPLE*10, t_end - t), d_in, &
                    GC_TOL, dt_gc, entered)
                dt = dt_gc
                mu_now = gc%mu
                if (entered) then
                    call gc_to_fo(pm, gc, lift, pm%p_phi, z, pphi, ok)
                    if (.not. ok) nfail = nfail + 1
                    pm%p_phi = pphi
                    region = region_of(pm, z(1:2))
                    e0 = hamiltonian(pm, region, z)
                    in_fo = .true.
                    nswitch = nswitch + 1
                end if
            end if
            t = t + dt
            if (t <= T_AVG) m0 = m0 + mu_now*dt
            if (t > t_end - T_AVG) m1 = m1 + mu_now*dt
        end do
        dp_phi = pm%p_phi - pphi0
        if (.not. in_fo) then
            call gc_to_fo(pm, gc, lift, pm%p_phi, z, pphi, ok)
            pm%p_phi = pphi
        end if
        region = region_of(pm, z(1:2))
        d_e = hamiltonian(pm, region, z)/hamiltonian(par_at(pphi0), 1, &
            initial_z(im)) - 1.0_dp
        frac = t_fo/t
    end subroutine run_marker

    real(dp) function sheet_s(gc) result(s)
        type(gc_state_t), intent(in) :: gc

        s = hypot(gc%r - par%r_major, gc%z) - par%r_sheet
    end function sheet_s

    function par_at(pphi) result(p)
        real(dp), intent(in) :: pphi
        type(field_params_t) :: p

        p = par
        p%p_phi = pphi
    end function par_at

    function initial_z(im) result(z)
        integer, intent(in) :: im
        real(dp) :: z(4)

        call initial_state(par, speed, mod(im - 1, 8) + 1, z)
        z(3) = z(3)*(1.0_dp + 1.0e-3_dp*real(im, dp)/NMARKER)
    end function initial_z

end program study_hybrid
