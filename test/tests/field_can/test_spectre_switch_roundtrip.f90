program test_spectre_switch_roundtrip
    ! Diagnostic: zeroth-order GC -> FO -> GC switch round trip with no FO step.
    ! spectre_fo_enter places the particle at x_gc + rho_L(alpha); the gyrophase
    ! sweep reproduces that construction for every alpha (alpha = 0 is checked
    ! against spectre_fo_enter itself). spectre_fo_to_gc then maps straight
    ! back. A second scan runs the full orbit between the two maps and samples
    ! the exit gyrophase, which is what a real excursion does. The defects
    ! d(mu)/mu, d(v_par)/v, d(x_gc) are recorded, not bounded: the test only
    ! requires the maps to complete with finite output.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use util, only: c, e_charge, p_mass, ev
    use magfie_sub, only: spectre_field, set_magfie_spectre_field
    use field_spectre, only: spectre_field_t, create_spectre_field
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_configure, &
        spectre_fo_enter, spectre_fo_to_gc, spectre_fo_radial_larmor, &
        spectre_fo_advance_until_exit, SPECTRE_FO_OK
    implicit none

    integer, parameter :: NPH = 16, NE = 3, NPOS = 5, NANG = 3, NLAM = 3
    real(dp), parameter :: PI = 3.14159265358979323846_dp
    ! Alpha particles at 3.5 MeV/facE_al with facE_al = 1000 (the NCSX study),
    ! and a factor 4 below and above in energy.
    real(dp), parameter :: E_SCALE(NE) = [0.25_dp, 1.0_dp, 4.0_dp]
    real(dp), parameter :: LAMBDAS(NLAM) = [-0.8_dp, 0.1_dp, 0.5_dp]
    ! Positions in radial gyroradii from interface 1: entry (3) and exit (4)
    ! distances of the default hybrid, both sides, and a mid-volume control.
    real(dp), parameter :: DIST(NPOS) = [-4.0_dp, -3.0_dp, 3.0_dp, 4.0_dp, 0.0_dp]
    character(len=10), parameter :: POS_NAME(NPOS) = [character(len=10) :: &
        'in -4 rL', 'in -3 rL', 'out +3 rL', 'out +4 rL', 'rho=0.5']
    type(spectre_field_t) :: sf
    type(spectre_fo_state_t) :: state, ref
    character(len=1024) :: h5file
    real(dp) :: ro0, ro0_bar, theta, zeta, lam, alpha, rho_r, rho_l
    real(dp) :: u(3), x_gc(3), b(3), bmod, e1(3), e2(3), eperp(3)
    real(dp) :: vpar, vperp, y(5), y_out(5), mu_in, mu_out, x_out(3), b_out(3)
    real(dp) :: bmod_out, xp(3), bp(3), bmod_p
    real(dp) :: dmu, dlam, drho, dang, dx, eps_loc
    real(dp) :: s_mu(NE, NPOS), s_mu1(NE, NPOS), s_lam(NE, NPOS)
    real(dp) :: s_rho(NE, NPOS), s_ang(NE, NPOS), s_x(NE, NPOS)
    real(dp) :: s_eps(NE, NPOS), max_mu(NE, NPOS), rl_cm(NE)
    integer :: n(NE, NPOS), ie, ip, ia, il, k, status, owner, owner_out, unit
    integer :: ncross(NE, NPOS), nfail(NE, NPOS)

    call get_command_argument(1, h5file)
    call create_spectre_field(sf, trim(h5file), status)
    if (status /= 0) error stop 'fixture initialization failed'
    call set_magfie_spectre_field(sf)
    if (spectre_field%data%Mvol < 2) error stop 'fixture needs an interface'
    call spectre_fo_configure(1.0e-9_dp)

    s_mu = 0; s_mu1 = 0; s_lam = 0; s_rho = 0; s_ang = 0; s_x = 0; s_eps = 0
    max_mu = 0; n = 0; ncross = 0; nfail = 0; rl_cm = 0
    open (newunit=unit, file='switch_roundtrip.dat', status='replace', &
          action='write')
    write (unit, '(A)') '# escale pos theta zeta lambda alpha owner_in '// &
        'owner_particle dmu_over_mu dlambda drho_over_rhor dtheta dzeta '// &
        'dxgc_over_rhoL eps_loc'
    do ie = 1, NE
        ro0 = alpha_ro0(3.5e6_dp/1000.0_dp*E_SCALE(ie))
        ro0_bar = ro0/sqrt(2.0_dp)
        do ip = 1, NPOS
            do ia = 1, NANG
                theta = real(ia - 1, dp)*PI/2.0_dp
                zeta = real(ia - 1, dp)*PI/9.0_dp
                do il = 1, NLAM
                    lam = LAMBDAS(il)
                    owner = merge(2, 1, DIST(ip) > 0.0_dp)
                    y = [merge(0.5_dp, 1.0_dp, ip == NPOS), theta, zeta, &
                         1.0_dp, lam]
                    call spectre_fo_radial_larmor(y, owner, ro0_bar, rho_r, &
                                                  status)
                    if (status /= SPECTRE_FO_OK) error stop 'larmor sample'
                    if (ip /= NPOS) then
                        y(1) = 1.0_dp + DIST(ip)*rho_r
                        call spectre_fo_radial_larmor(y, owner, ro0_bar, &
                                                      rho_r, status)
                        if (status /= SPECTRE_FO_OK) error stop 'larmor sample'
                        y(1) = 1.0_dp + DIST(ip)*rho_r
                    end if
                    u = y(1:3)
                    call sample(u, x_gc, b, bmod)
                    vpar = sqrt(2.0_dp)*lam
                    vperp = sqrt(2.0_dp)*sqrt(1.0_dp - lam**2)
                    rho_l = ro0_bar*vperp/bmod
                    rl_cm(ie) = max(rl_cm(ie), rho_l)
                    mu_in = vperp**2/(2.0_dp*bmod)
                    call perpendicular_frame(b, e1, e2)
                    call spectre_fo_enter(ref, y, owner, 1, ro0_bar, status)
                    if (status /= SPECTRE_FO_OK) error stop 'spectre_fo_enter'
                    do k = 1, NPH
                        alpha = 2.0_dp*PI*real(k - 1, dp)/real(NPH, dp)
                        eperp = cos(alpha)*e1 + sin(alpha)*e2
                        state = spectre_fo_state_t()
                        state%v = vpar*b + vperp*eperp
                        state%x = x_gc + (ro0_bar/bmod)*cross(b, vperp*eperp)
                        state%u = u
                        state%owner = owner
                        state%entry_owner = owner
                        state%iface = 1
                        state%active = .true.
                        if (k == 1) then
                            if (norm2(state%x - ref%x) > 1.0e-9_dp*rho_l .or. &
                                norm2(state%v - ref%v) > 1.0e-12_dp) &
                                error stop 'phase sweep differs from enter'
                            state = ref
                        end if
                        call spectre_fo_to_gc(state, ro0_bar, y_out, &
                                              owner_out, status, mu_out)
                        if (status /= SPECTRE_FO_OK) then
                            ! A failed production map is itself recorded.
                            nfail(ie, ip) = nfail(ie, ip) + 1
                            cycle
                        end if
                        call sample(state%u, xp, bp, bmod_p)
                        call sample(y_out(1:3), x_out, b_out, bmod_out)
                        dmu = mu_out/mu_in - 1.0_dp
                        dlam = y_out(5) - lam
                        drho = (y_out(1) - y(1))/rho_r
                        dang = sqrt(wrap(y_out(2) - y(2))**2 + &
                            wrap(y_out(3) - y(3))**2)
                        dx = norm2(x_out - x_gc)/rho_l
                        eps_loc = abs(bmod_p - bmod)/bmod
                        if (.not. all(ieee_is_finite([dmu, dlam, drho, dang, &
                            dx, y_out(4)]))) error stop 'non-finite defect'
                        if (abs(y_out(4) - 1.0_dp) > 1.0e-12_dp) &
                            error stop 'switch changed the speed'
                        write (unit, '(F6.2,I3,4F9.4,2I3,7ES13.4)') E_SCALE(ie), &
                            ip, theta, zeta, lam, alpha, owner, state%owner, &
                            dmu, dlam, drho, wrap(y_out(2) - y(2)), &
                            wrap(y_out(3) - y(3)), &
                            dx, eps_loc
                        if (state%owner /= owner) &
                            ncross(ie, ip) = ncross(ie, ip) + 1
                        n(ie, ip) = n(ie, ip) + 1
                        s_mu(ie, ip) = s_mu(ie, ip) + dmu**2
                        s_mu1(ie, ip) = s_mu1(ie, ip) + dmu
                        s_lam(ie, ip) = s_lam(ie, ip) + dlam**2
                        s_rho(ie, ip) = s_rho(ie, ip) + drho**2
                        s_ang(ie, ip) = s_ang(ie, ip) + dang**2
                        s_x(ie, ip) = s_x(ie, ip) + dx**2
                        s_eps(ie, ip) = s_eps(ie, ip) + eps_loc**2
                        max_mu(ie, ip) = max(max_mu(ie, ip), abs(dmu))
                    end do
                end do
            end do
        end do
    end do
    close (unit)

    print '(A)', 'GC->FO->GC round trip, no FO step; rms over theta, zeta, '// &
        'lambda and 16 gyrophases'
    print '(A)', ' E/E0  rhoL[cm] position    n ncross nfail rms(dmu/mu) '// &
        'mean(dmu/mu) max|dmu/mu|  rms(dlam)  rms(drho/rhor) rms(dang) '// &
        'rms(dxgc/rhoL) rms(eps_loc)'
    do ip = 1, NPOS
        do ie = 1, NE
            print '(F5.2,F9.4,1X,A10,3I6,8ES12.3)', E_SCALE(ie), rl_cm(ie), &
                POS_NAME(ip), n(ie, ip), ncross(ie, ip), nfail(ie, ip), &
                sqrt(s_mu(ie, ip)/n(ie, ip)), s_mu1(ie, ip)/n(ie, ip), &
                max_mu(ie, ip), sqrt(s_lam(ie, ip)/n(ie, ip)), &
                sqrt(s_rho(ie, ip)/n(ie, ip)), sqrt(s_ang(ie, ip)/n(ie, ip)), &
                sqrt(s_x(ie, ip)/n(ie, ip)), sqrt(s_eps(ie, ip)/n(ie, ip))
        end do
    end do
    print '(A)', 'scaling exponent of rms(dmu/mu) and |mean(dmu/mu)| in '// &
        'gyroradius (E x4 = rhoL x2):'
    do ip = 1, NPOS
        print '(1X,A10,A,2F7.3,A,2F7.3)', POS_NAME(ip), '  rms: ', &
            (log(sqrt(s_mu(ie + 1, ip)/s_mu(ie, ip)))/log(2.0_dp), ie=1, NE - 1), &
            '  mean: ', (log(abs(s_mu1(ie + 1, ip)/s_mu1(ie, ip)))/log(2.0_dp), &
                         ie=1, NE - 1)
    end do
    if (any(n + nfail /= NANG*NLAM*NPH)) error stop 'missing round trips'
    if (any(n == 0)) error stop 'no completed round trip in a class'
    if (.not. all(ieee_is_finite(s_mu))) error stop 'non-finite summary'

    call exit_phase_scan

contains

    subroutine exit_phase_scan
        ! Enter at alpha = 0, then run the full orbit (smooth field, no sheet
        ! reached) and apply spectre_fo_to_gc at 16 exit gyrophases per
        ! gyroperiod over two periods. The spread of the returned mu over exit
        ! phase is the per-excursion switch defect of the zeroth-order pair;
        ! the no-step round trip above cancels it by symmetry.
        integer, parameter :: NSAMP = 32
        real(dp) :: mus(NSAMP), lams(NSAMP), rhos(NSAMP), t_gyro, used, yy(5)
        real(dp) :: amp_mu(NE, 2), rms_mu(NE, 2), mean_mu(NE, 2), msq_mu(NE, 2)
        real(dp) :: amp_lam(NE, 2), amp_rho(NE, 2), eps_r(NE, 2)
        integer :: nc(NE, 2), ic, j, own, nbad(NE, 2)
        logical :: exited

        amp_mu = 0; rms_mu = 0; mean_mu = 0; msq_mu = 0; amp_lam = 0; amp_rho = 0
        eps_r = 0; nc = 0; nbad = 0
        do ie = 1, NE
            ro0 = alpha_ro0(3.5e6_dp/1000.0_dp*E_SCALE(ie))
            ro0_bar = ro0/sqrt(2.0_dp)
            do ic = 1, 2
                do ia = 1, NANG
                    do il = 1, NLAM
                        lam = LAMBDAS(il)
                        y = [0.5_dp, real(ia - 1, dp)*PI/2.0_dp, &
                             real(ia - 1, dp)*PI/9.0_dp, 1.0_dp, lam]
                        if (ic == 1) then
                            call spectre_fo_radial_larmor(y, 1, ro0_bar, rho_r, &
                                                          status)
                            y(1) = 1.0_dp - 3.0_dp*rho_r
                        end if
                        call spectre_fo_radial_larmor(y, 1, ro0_bar, rho_r, status)
                        call sample(y(1:3), x_gc, b, bmod)
                        mu_in = (1.0_dp - lam**2)/bmod
                        call spectre_fo_enter(state, y, 1, 1, ro0_bar, status)
                        if (status /= SPECTRE_FO_OK) error stop 'scan enter'
                        state%exit_rho = 10.0_dp
                        t_gyro = 2.0_dp*PI*ro0_bar/bmod
                        do j = 1, NSAMP
                            call spectre_fo_advance_until_exit(state, &
                                t_gyro/16.0_dp, ro0_bar, used, yy, own, &
                                exited, status)
                            if (status /= SPECTRE_FO_OK .or. exited) &
                                error stop 'scan full orbit failed'
                            call spectre_fo_to_gc(state, ro0_bar, yy, own, &
                                                  status, mu_out)
                            if (status /= SPECTRE_FO_OK) exit
                            mus(j) = mu_out/mu_in - 1.0_dp
                            lams(j) = yy(5) - lam
                            rhos(j) = (yy(1) - y(1))/rho_r
                        end do
                        if (status /= SPECTRE_FO_OK) then
                            nbad(ie, ic) = nbad(ie, ic) + 1
                            cycle
                        end if
                        if (.not. all(ieee_is_finite(mus))) &
                            error stop 'non-finite scan defect'
                        nc(ie, ic) = nc(ie, ic) + 1
                        amp_mu(ie, ic) = amp_mu(ie, ic) + &
                                         (0.5_dp*(maxval(mus) - minval(mus)))**2
                        rms_mu(ie, ic) = rms_mu(ie, ic) + sum(mus**2)/NSAMP
                        mean_mu(ie, ic) = mean_mu(ie, ic) + sum(mus)/NSAMP
                        msq_mu(ie, ic) = msq_mu(ie, ic) + (sum(mus)/NSAMP)**2
                        amp_lam(ie, ic) = amp_lam(ie, ic) + &
                                          (0.5_dp*(maxval(lams) - minval(lams)))**2
                        amp_rho(ie, ic) = amp_rho(ie, ic) + &
                                          (0.5_dp*(maxval(rhos) - minval(rhos)))**2
                        eps_r(ie, ic) = eps_r(ie, ic) + &
                                        (ro0_bar*sqrt(2.0_dp*(1.0_dp - lam**2))/bmod)**2
                    end do
                end do
            end do
        end do
        print '(A)', 'Exit-phase scan: enter at alpha=0, full orbit, to_gc at '// &
            '32 exit phases over 2 gyroperiods (rms over theta, zeta, lambda)'
        print '(A)', ' E/E0 position    n nfail  amp(dmu/mu) rms(dmu/mu) '// &
            'mean(dmu/mu) rms(mean)   amp(dlam)  amp(drho/rhor) rhoL/R0'
        do ic = 1, 2
            do ie = 1, NE
                if (nc(ie, ic) == 0) error stop 'no completed exit-phase scan'
                print '(F5.2,1X,A10,2I5,7ES12.3)', E_SCALE(ie), &
                    merge('in -3 rL  ', 'rho=0.5   ', ic == 1), nc(ie, ic), &
                    nbad(ie, ic), sqrt(amp_mu(ie, ic)/nc(ie, ic)), &
                    sqrt(rms_mu(ie, ic)/nc(ie, ic)), mean_mu(ie, ic)/nc(ie, ic), &
                    sqrt(msq_mu(ie, ic)/nc(ie, ic)), sqrt(amp_lam(ie, ic)/nc(ie, ic)), &
                    sqrt(amp_rho(ie, ic)/nc(ie, ic)), &
                    sqrt(eps_r(ie, ic)/nc(ie, ic))/139.5_dp
            end do
            print '(A,2F7.3)', '   scaling exponent of amp(dmu/mu): ', &
                (log(sqrt(amp_mu(ie + 1, ic)/amp_mu(ie, ic)))/log(2.0_dp), &
                 ie=1, NE - 1)
        end do
    end subroutine exit_phase_scan

    pure real(dp) function wrap(d) result(w)
        real(dp), intent(in) :: d
        w = d - 2.0_dp*PI*real(nint(d/(2.0_dp*PI)), dp)
    end function wrap

    real(dp) function alpha_ro0(energy_ev) result(ro)
        !> SIMPLE's rlarm = v0 m c / (Z e) for an alpha particle (params.f90).
        real(dp), intent(in) :: energy_ev
        real(dp) :: v0
        v0 = sqrt(2.0_dp*energy_ev*ev/(4.0_dp*p_mass))
        ro = v0*4.0_dp*p_mass*c/(2.0_dp*e_charge)
    end function alpha_ro0

    subroutine sample(u, x, b, bmod)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3), b(3), bmod
        real(dp) :: a(3), hc(3), sqgb(3), g(3, 3), gi(3, 3), sg, e(3, 3)

        call spectre_field%coords%evaluate_cart(u, x)
        call spectre_field%coords%covariant_basis(u, e)
        call spectre_field%coords%metric_tensor(u, g, gi, sg)
        call spectre_field%evaluate(u, a, hc, bmod, sqgb)
        b = matmul(e, sqgb/sg)/bmod
        x = x*1.0e2_dp
        bmod = bmod*1.0e4_dp
    end subroutine sample

    pure subroutine perpendicular_frame(b, e1, e2)
        real(dp), intent(in) :: b(3)
        real(dp), intent(out) :: e1(3), e2(3)
        real(dp) :: axis(3)

        if (abs(b(3)) < 0.9_dp) then
            axis = [0.0_dp, 0.0_dp, 1.0_dp]
        else
            axis = [1.0_dp, 0.0_dp, 0.0_dp]
        end if
        e1 = axis - dot_product(axis, b)*b
        e1 = e1/norm2(e1)
        e2 = cross(b, e1)
    end subroutine perpendicular_frame

    pure function cross(a, b) result(c3)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c3(3)
        c3 = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross

end program test_spectre_switch_roundtrip
