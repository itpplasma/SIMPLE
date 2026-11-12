program test_spectre_switch_mu_bias
    ! Guiding centre -> full orbit -> guiding centre passages near the
    ! current-carrying interface of the axisymmetric tok2vol fixture. The
    ! instantaneous first-order switches change the guiding-centre mu by an
    ! O(rho*) amount whose mean over gyrophase does not vanish at a sheet. The
    ! switch error of a passage is the guiding-centre change of mu minus the
    ! change of the full orbit's gyroaveraged mu over it; its mean must stay
    ! below MU_BIAS_BOUND, and the guiding-centre p_zeta of entry and exit
    ! must agree.
    ! With the argument "compare" the instantaneous switches and the averaged
    ! exit without the p_zeta shift are measured as well (no bounds applied).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_sub, only: spectre_field, set_magfie_spectre_field
    use field_spectre, only: spectre_field_t, create_spectre_field
    use spectre_fo_field, only: spectre_fo_field_local, spectre_fo_local_t, &
                                spectre_fo_field_locate
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_configure, &
                                 spectre_fo_enter, spectre_fo_advance_until_exit, &
                                 spectre_fo_radial_larmor, SPECTRE_FO_OK
    implicit none

    ! Radial gyroradius at the interface in units of the volume width.
    real(dp), parameter :: RHO_STAR = 0.02_dp
    ! Entry and exit distances of the default hybrid (spectre_fo_entry_larmor).
    real(dp), parameter :: ENTRY_LARMOR = 1.5_dp
    integer, parameter :: NTH = 16, NLAM = 4
    real(dp), parameter :: LAMBDAS(NLAM) = [-0.6_dp, -0.2_dp, 0.2_dp, 0.6_dp]
    ! Measured with "compare": instantaneous switches -3.2e-4 (rms 9.4e-4),
    ! averaged -1.4e-11 (rms 1.2e-10).
    real(dp), parameter :: MU_BIAS_BOUND = 1.0e-5_dp
    integer, parameter :: MAX_GYRO = 2000
    real(dp), parameter :: PZETA_TOL = 1.0e-10_dp
    real(dp), parameter :: PI = 3.14159265358979323846_dp

    type(spectre_field_t) :: sf
    character(len=1024) :: h5file, mode
    real(dp) :: ro0_bar, mean, rms, dpz
    integer :: status, n, nfb

    call get_command_argument(1, h5file)
    mode = ''
    if (command_argument_count() > 1) call get_command_argument(2, mode)
    call create_spectre_field(sf, trim(h5file), status)
    if (status /= 0) error stop 'fixture initialization failed'
    call set_magfie_spectre_field(sf)
    if (spectre_field%data%Mvol < 2) error stop 'fixture needs an interface'
    call spectre_fo_configure(1.0e-9_dp)
    ro0_bar = scaled_ro0()

    if (trim(mode) == 'compare') then
        call spectre_fo_configure(1.0e-9_dp, averaged=.false.)
        call passages(ro0_bar, mean, rms, dpz, n, nfb)
        call report('instantaneous', mean, rms, dpz, n, nfb)
        call spectre_fo_configure(1.0e-9_dp, pzeta_shift=.false.)
        call passages(ro0_bar, mean, rms, dpz, n, nfb)
        call report('averaged, no p_zeta shift', mean, rms, dpz, n, nfb)
        call spectre_fo_configure(1.0e-9_dp)
    end if
    call passages(ro0_bar, mean, rms, dpz, n, nfb)
    call report('averaged', mean, rms, dpz, n, nfb)
    if (n < NTH*NLAM/4) error stop 'too few completed passages'
    if (abs(mean) > MU_BIAS_BOUND) error stop 'mean mu change per passage too large'
    if (dpz > PZETA_TOL) error stop 'p_zeta not conserved over the passage'

contains

    real(dp) function scaled_ro0() result(ro0_bar)
        ! Choose ro0 such that a lambda = 0 guiding centre on the interface
        ! has the radial gyroradius RHO_STAR (it is linear in ro0).
        real(dp) :: rho_r, worst

        worst = 0.0_dp
        do n = 1, NTH
            call spectre_fo_radial_larmor([1.0_dp, 2.0_dp*PI*(n - 1)/NTH, &
                                      0.0_dp, 1.0_dp, 0.0_dp], 1, 1.0_dp, rho_r, status)
            if (status /= SPECTRE_FO_OK) error stop 'larmor sample'
            worst = max(worst, rho_r)
        end do
        ro0_bar = RHO_STAR/worst
    end function scaled_ro0

    subroutine passages(ro0_bar, mean, rms, dpz_max, n, nfallback)
        !> Switch error per passage: the guiding-centre change of mu minus the
        !> change of the full orbit's gyroaveraged mu over the same passage,
        !> each average taken over the gyroperiod after entry and after exit.
        real(dp), intent(in) :: ro0_bar
        real(dp), intent(out) :: mean, rms, dpz_max
        integer, intent(out) :: n, nfallback
        type(spectre_fo_state_t) :: state
        real(dp) :: y(5), y_out(5), rho_r, used, err, t_chunk, gauge
        real(dp) :: p_in, p_out, mu_fo_in, mu_fo_out, raw
        integer :: it, il, owner, st, ic
        logical :: exited

        gauge = acov_zeta(2) - acov_zeta(1)
        mean = 0.0_dp
        rms = 0.0_dp
        raw = 0.0_dp
        dpz_max = 0.0_dp
        n = 0
        nfallback = 0
        do it = 1, NTH
            do il = 1, NLAM
                y = [1.0_dp, 2.0_dp*PI*(it - 1)/NTH, 0.0_dp, 1.0_dp, LAMBDAS(il)]
                call spectre_fo_radial_larmor(y, 1, ro0_bar, rho_r, st)
                y(1) = 1.0_dp - ENTRY_LARMOR*rho_r
                call spectre_fo_enter(state, y, 1, 1, ro0_bar, st)
                if (st /= SPECTRE_FO_OK) error stop 'entry failed'
                if (state%entry_fallback) nfallback = nfallback + 1
                state%exit_rho = (ENTRY_LARMOR + 1.0_dp)*rho_r
                mu_fo_in = averaged_mu(state, ro0_bar)
                t_chunk = 100.0_dp*PI*ro0_bar/b_local(y, 1)
                do ic = 1, MAX_GYRO/50
                    call spectre_fo_advance_until_exit(state, t_chunk, ro0_bar, &
                                                       used, y_out, owner, exited, st)
                    if (st /= SPECTRE_FO_OK) error stop 'full orbit failed'
                    if (exited) exit
                end do
                if (.not. exited) cycle
                if (state%exit_fallback) nfallback = nfallback + 1
                mu_fo_out = averaged_mu(state, ro0_bar)
                err = gc_mu(y_out, owner)/gc_mu(y, 1) - mu_fo_out/mu_fo_in
                p_in = gc_pzeta(y, 1, ro0_bar)
                p_out = gc_pzeta(y_out, owner, ro0_bar) &
                        - merge(gauge/ro0_bar, 0.0_dp, owner == 2)
                dpz_max = max(dpz_max, abs(p_out - p_in)/max(abs(p_in), 1.0_dp))
                n = n + 1
                mean = mean + err
                rms = rms + err**2
                raw = raw + gc_mu(y_out, owner)/gc_mu(y, 1) - 1.0_dp
            end do
        end do
        mean = mean/max(n, 1)
        rms = sqrt(rms/max(n, 1))
        print '(A,ES11.3)', 'mean guiding-centre dmu/mu per passage (switches '// &
            'and physics): ', raw/max(n, 1)
    end subroutine passages

    real(dp) function averaged_mu(start, ro0_bar) result(mu)
        !> mu = <v_perp^2>/(2 B(<X>)) over one gyroperiod 2 pi ro0/B(X) of a
        !> copy of the full orbit, with the first-order X = x - (ro0/B) b x v,
        !> v_perp^2 = v^2 - (v.b(X))^2, at NAVG equally spaced times.
        type(spectre_fo_state_t), intent(in) :: start
        real(dp), intent(in) :: ro0_bar
        integer, parameter :: NAVG = 32
        type(spectre_fo_state_t) :: s
        type(spectre_fo_local_t) :: gc
        real(dp) :: x_sum(3), vperp2, dt, used, y(5), u(3), xg(3)
        integer :: k, owner, st
        logical :: exited

        s = start
        s%active = .true.
        s%exit_rho = huge(1.0_dp)
        call first_order(s, ro0_bar, u, owner, xg, gc)
        dt = 2.0_dp*PI*ro0_bar/gc%bmod/NAVG
        x_sum = 0.0_dp
        vperp2 = 0.0_dp
        do k = 1, NAVG
            call spectre_fo_advance_until_exit(s, dt, ro0_bar, used, y, owner, &
                                               exited, st)
            if (st /= SPECTRE_FO_OK .or. exited) error stop 'averaging orbit failed'
            call first_order(s, ro0_bar, u, owner, xg, gc)
            x_sum = x_sum + xg
            vperp2 = vperp2 + dot_product(s%v, s%v) - dot_product(s%v, gc%b)**2
        end do
        call spectre_fo_field_locate(x_sum/NAVG, u, owner, st)
        if (st /= 0) error stop 'averaged centre not located'
        call spectre_fo_field_local(u, owner, gc)
        mu = vperp2/NAVG/(2.0_dp*gc%bmod)
    end function averaged_mu

    subroutine first_order(s, ro0_bar, u, owner, xg, gc)
        type(spectre_fo_state_t), intent(in) :: s
        real(dp), intent(in) :: ro0_bar
        real(dp), intent(out) :: u(3), xg(3)
        integer, intent(out) :: owner
        type(spectre_fo_local_t), intent(out) :: gc
        type(spectre_fo_local_t) :: p
        integer :: st

        call spectre_fo_field_local(s%u, s%owner, p)
        xg = p%x - (ro0_bar/p%bmod)*cross(p%b, s%v)
        u = s%u
        owner = s%owner
        call spectre_fo_field_locate(xg, u, owner, st)
        if (st /= 0) error stop 'guiding centre not located'
        call spectre_fo_field_local(u, owner, gc)
    end subroutine first_order

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross

    subroutine report(label, mean, rms, dpz, n, nfb)
        character(*), intent(in) :: label
        real(dp), intent(in) :: mean, rms, dpz
        integer, intent(in) :: n, nfb
        print '(A,T30,A,I4,A,ES11.3,A,ES11.3,A,ES11.3,A,I3)', label, &
            ' passages', n, '  switch dmu/mu mean', mean, '  rms', rms, &
            '  max dp_zeta', dpz, '  instantaneous', nfb
    end subroutine report

    real(dp) function gc_mu(y, owner) result(mu)
        real(dp), intent(in) :: y(5)
        integer, intent(in) :: owner
        type(spectre_fo_local_t) :: loc
        call spectre_fo_field_local(y(1:3), owner, loc)
        mu = y(4)**2*(1.0_dp - y(5)**2)/loc%bmod
    end function gc_mu

    real(dp) function gc_pzeta(y, owner, ro0_bar) result(p)
        real(dp), intent(in) :: y(5), ro0_bar
        integer, intent(in) :: owner
        type(spectre_fo_local_t) :: loc
        call spectre_fo_field_local(y(1:3), owner, loc)
        p = sqrt(2.0_dp)*y(4)*y(5)*dot_product(loc%b, loc%e(:, 3)) &
            + loc%acov(3)/ro0_bar
    end function gc_pzeta

    real(dp) function b_local(y, owner) result(bmod)
        real(dp), intent(in) :: y(5)
        integer, intent(in) :: owner
        type(spectre_fo_local_t) :: loc
        call spectre_fo_field_local(y(1:3), owner, loc)
        bmod = loc%bmod
    end function b_local

    real(dp) function acov_zeta(owner) result(a)
        ! A_zeta on the interface, constant in axisymmetry, in volume owner's
        ! gauge.
        integer, intent(in) :: owner
        type(spectre_fo_local_t) :: loc
        call spectre_fo_field_local([1.0_dp, 0.0_dp, 0.0_dp], owner, loc)
        a = loc%acov(3)
    end function acov_zeta

end program test_spectre_switch_mu_bias
