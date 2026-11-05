module collis_alp
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none

    integer, parameter :: nsorts = 3
    real(wp) :: collision_s_max = 1.0_wp

    real(wp), dimension(nsorts) :: efcolf, velrat, enrat
    ! Coulomb logarithm at the reference speed and the ratio of classical to
    ! quantum minimum impact parameter there; a zero ratio keeps it constant.
    real(wp), dimension(nsorts) :: clog_ref = 0.0_wp, clog_bratio = 0.0_wp
    ! Species of the radial collision profile, set by init_collision_profiles.
    real(wp) :: profile_am1 = 0.0_wp, profile_am2 = 0.0_wp
    real(wp) :: profile_Z1 = 0.0_wp, profile_Z2 = 0.0_wp
    real(wp) :: profile_ealpha = 0.0_wp
    real(wp) :: profile_test_mass = 4.0_wp, profile_test_charge = 2.0_wp
    logical :: profile_ready = .false.
    ! Radially uniform coefficients override the profile when set.
    logical :: uniform_ready = .false.
    real(wp), dimension(nsorts) :: uniform_efcolf = 0.0_wp, &
                                   uniform_velrat = 0.0_wp, uniform_enrat = 0.0_wp, &
                                   uniform_clog = 0.0_wp, uniform_bratio = 0.0_wp
    real(wp) :: lorentz_nu = 0.0_wp
    real(wp) :: lorentz_dhh = 0.0_wp

contains

    subroutine configure_lorentz_nustar(nu_star, major_radius_cm, iota, vref)
        real(wp), intent(in) :: nu_star, major_radius_cm, iota, vref

        if (nu_star < 0.0_wp) error stop 'nu_star_standard must be nonnegative'
        if (major_radius_cm <= 0.0_wp) &
            error stop 'lorentz_major_radius_cm must be positive'
        if (abs(iota) <= tiny(1.0_wp)) error stop 'lorentz_iota must be nonzero'
        if (vref <= 0.0_wp) error stop 'Lorentz reference velocity must be positive'

        lorentz_nu = nu_star*abs(iota)*vref/major_radius_cm
        ! stost advances in path length d tau = vref dt.  This coefficient gives
        ! the Lorentz generator (nu/2) d/dlambda[(1-lambda^2)d/dlambda].
        lorentz_dhh = lorentz_nu/(2.0_wp*vref)
    end subroutine configure_lorentz_nustar

    subroutine stost_lorentz(z, dtauc, ierr)
        real(wp), intent(inout) :: z(5)
        real(wp), intent(in) :: dtauc
        integer, intent(out) :: ierr

        real(wp) :: alam, coala, dalam
        real :: ur

        ierr = 0
        alam = z(5)
        coala = 1.0_wp - alam**2
        if (coala < 0.0_wp) then
            ierr = 1
            return
        end if

        call getran(1, ur)
        dalam = sqrt(2.0_wp*lorentz_dhh*coala*dtauc)*real(ur, wp) &
                - 2.0_wp*alam*lorentz_dhh*dtauc
        alam = alam + dalam
        if (alam > 1.0_wp) then
            alam = 2.0_wp - alam
        else if (alam < -1.0_wp) then
            alam = -2.0_wp - alam
        end if
        if (abs(alam) > 1.0_wp) then
            ierr = 2
            return
        end if
        z(5) = alam
    end subroutine stost_lorentz

    subroutine coleff(p, dpp, dhh, fpeff)
        real(wp), intent(in) :: p
        real(wp), intent(out) :: dpp, dhh, fpeff

        call coleff_local(p, efcolf, velrat, enrat, dpp, dhh, fpeff, &
                          clog_ref, clog_bratio)
    end subroutine coleff

    subroutine coleff_local(p, efcolf_loc, velrat_loc, enrat_loc, dpp, dhh, fpeff, &
                            clog_loc, bratio_loc)
        real(wp), intent(in) :: p
       real(wp), intent(in) :: efcolf_loc(nsorts), velrat_loc(nsorts), enrat_loc(nsorts)
        real(wp), intent(out) :: dpp, dhh, fpeff
        real(wp), intent(in), optional :: clog_loc(nsorts), bratio_loc(nsorts)

        integer :: i
        real(wp) :: plim, xbeta, dp, dh, dpd, f, df

        plim = max(p, 1.d-8)

        dpp = 0.0d0
        dhh = 0.0d0
        fpeff = 0.0d0

        do i = 1, nsorts
            xbeta = p*velrat_loc(i)
            call onseff(xbeta, dp, dh, dpd)
            call log_factor(i, p, velrat_loc, f, df, clog_loc, bratio_loc)
            dpp = dpp + f*dp*efcolf_loc(i)
            dhh = dhh + f*dh*efcolf_loc(i)
            ! The speed derivative of the logarithm enters the Ito drift.
            fpeff = fpeff + (f*(dpd/plim - 2.0d0*dp*p*enrat_loc(i)) + df*dp) &
                    *efcolf_loc(i)
        end do

        dhh = dhh/plim**2
    end subroutine coleff_local

    subroutine onseff(v, dp, dh, dpd)
        real(wp), intent(in) :: v
        real(wp), intent(out) :: dp, dh, dpd

        real(wp), parameter :: sqp = 1.7724538d0
        real(wp), parameter :: cons = 0.75225278d0
        real(wp) :: v2, v3, ex, er

        v2 = v**2
        v3 = v2*v
        if (v < 0.01d0) then
            dp = cons*(1.d0 - 0.6d0*v2)
            dh = cons*(1.d0 - 0.2d0*v2)
            dpd = 2.d0*cons*(1.d0 - 1.2d0*v2)
        else if (v > 6.d0) then
            dp = 1.d0/v3
            dh = (1.d0 - 0.5d0/v2)/v
            dpd = -1.d0/v3
        else
            ex = exp(-v2)/sqp
            er = erf(v)
            dp = er/v3 - 2.d0*ex/v2
            dh = er*(1.d0 - 0.5d0/v2)/v + ex/v2
            dpd = 4.d0*ex - dp
        end if
    end subroutine onseff

    !> Cartesian-velocity coefficients of the same collision operator, regular
    !> at p = 0: parallel diffusion dpar, perpendicular diffusion dperp and the
    !> drift a along the velocity, with fpeff = a + 2*dperp/p.
    subroutine coleff_cartesian(p, efcolf_loc, velrat_loc, enrat_loc, &
                                dpar, dperp, a, clog_loc, bratio_loc)
        real(wp), intent(in) :: p
       real(wp), intent(in) :: efcolf_loc(nsorts), velrat_loc(nsorts), enrat_loc(nsorts)
        real(wp), intent(out) :: dpar, dperp, a
        real(wp), intent(in), optional :: clog_loc(nsorts), bratio_loc(nsorts)

        real(wp), parameter :: cons = 0.75225278d0
        integer :: i
        real(wp) :: x, dp, dh, dpd, g, f, df

        dpar = 0.0_wp
        dperp = 0.0_wp
        a = 0.0_wp
        do i = 1, nsorts
            x = p*velrat_loc(i)
            call onseff(x, dp, dh, dpd)
            if (x < 0.01d0) then
                g = -2.0d0*cons*x
            else if (x > 6.0d0) then
                g = -2.0d0/x**2
            else
                g = (dpd - 2.0d0*dh)/x
            end if
            call log_factor(i, p, velrat_loc, f, df, clog_loc, bratio_loc)
            dpar = dpar + f*dp*efcolf_loc(i)
            dperp = dperp + f*dh*efcolf_loc(i)
            a = a + efcolf_loc(i)*(f*(velrat_loc(i)*g - 2.0d0*dp*p*enrat_loc(i)) &
                                   + df*dp)
        end do
    end subroutine coleff_cartesian

    !> Coulomb logarithm at speed p relative to the reference speed, f = L(p)/L(1),
    !> and df = (dL/dp)/L(1). L = ln(lambda_D/b_min) with b_min the larger of the
    !> classical and quantum impact parameters at the mean relative speed
    !> vbar**2 = v0**2*(p**2 + 1/velrat**2).
    subroutine log_factor(i, p, velrat_loc, f, df, clog_loc, bratio_loc)
        integer, intent(in) :: i
        real(wp), intent(in) :: p, velrat_loc(nsorts)
        real(wp), intent(out) :: f, df
        real(wp), intent(in), optional :: clog_loc(nsorts), bratio_loc(nsorts)

        real(wp) :: w, q, lmin, lmin1

        f = 1.0_wp
        df = 0.0_wp
        if (.not. (present(clog_loc) .and. present(bratio_loc))) return
        if (bratio_loc(i) <= 0.0_wp .or. velrat_loc(i) <= 0.0_wp) return
        w = p**2 + 1.0_wp/velrat_loc(i)**2
        q = w/(1.0_wp + 1.0_wp/velrat_loc(i)**2)
        ! Logarithms of b_min relative to the quantum value at p = 1.
        lmin1 = max(log(bratio_loc(i)), 0.0_wp)
        if (log(bratio_loc(i)) - log(q) >= -0.5_wp*log(q)) then
            lmin = log(bratio_loc(i)) - log(q)
            df = 2.0_wp*p/w
        else
            lmin = -0.5_wp*log(q)
            df = p/w
        end if
        f = (clog_loc(i) - (lmin - lmin1))/clog_loc(i)
        df = df/clog_loc(i)
    end subroutine log_factor

    subroutine loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, &
                            tempe, ealpha, v0, dchichi, slowrate, &
                            dchichi_norm, slowrate_norm, test_mass, test_charge)
        real(wp), intent(in) :: am1, am2, Z1, Z2, densi1, densi2
        real(wp), intent(in) :: tempi1, tempi2, tempe, ealpha
        real(wp), intent(in), optional :: test_mass, test_charge
        real(wp), intent(out) :: v0, dchichi, slowrate, dchichi_norm, slowrate_norm

        call compute_collision_coeffs(am1, am2, Z1, Z2, densi1, densi2, &
                                      tempi1, tempi2, tempe, ealpha, v0, &
                                      efcolf, velrat, enrat, test_mass, test_charge, &
                                      clog_ref, clog_bratio)

        call compute_rates(v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
    end subroutine loacol_alpha

    subroutine compute_collision_coeffs(am1, am2, Z1, Z2, densi1, densi2, &
                                        tempi1, tempi2, tempe, ealpha, v0, &
                                        efcolf_out, velrat_out, enrat_out, &
                                        test_mass, test_charge, clog_out, bratio_out)
        use util, only: pi, hbar, pmass => p_mass, emass => e_mass, &
                        e_cgs => e_charge, ev_cgs => ev
        real(wp), intent(in) :: am1, am2, Z1, Z2, densi1, densi2
        real(wp), intent(in) :: tempi1, tempi2, tempe, ealpha
        real(wp), intent(in), optional :: test_mass, test_charge
        real(wp), intent(out) :: v0
      real(wp), intent(out) :: efcolf_out(nsorts), velrat_out(nsorts), enrat_out(nsorts)
        real(wp), intent(out), optional :: clog_out(nsorts), bratio_out(nsorts)

        real(wp) :: dense, vti1, vti2, vte
        real(wp) :: alami1, alami2, alame, frecol_base, mass, charge
        real(wp) :: debye, field_mass(nsorts), field_charge(nsorts), vth2(nsorts)
        real(wp) :: clog(nsorts), bratio(nsorts), reduced, vbar2, bcl, bqm
        integer :: i

        mass = 4.0_wp
        charge = 2.0_wp
        if (present(test_mass)) mass = test_mass
        if (present(test_charge)) charge = test_charge
        if (mass <= 0.0_wp .or. charge <= 0.0_wp) &
            error stop "Collision test species must have positive mass and charge"

        enrat_out(1) = ealpha/tempi1
        enrat_out(2) = ealpha/tempi2
        enrat_out(3) = ealpha/tempe

        v0 = sqrt(2.d0*ealpha*ev_cgs/(mass*pmass))
        vti1 = sqrt(2.d0*tempi1*ev_cgs/(pmass*am1))
        vti2 = sqrt(2.d0*tempi2*ev_cgs/(pmass*am2))
        vte = sqrt(2.d0*tempe*ev_cgs/emass)

        velrat_out(1) = v0/vti1
        velrat_out(2) = v0/vti2
        velrat_out(3) = v0/vte

        dense = densi1*Z1 + densi2*Z2
        if (dense == 0.0_wp) then
            efcolf_out = 0.0_wp
            if (present(clog_out)) clog_out = 0.0_wp
            if (present(bratio_out)) bratio_out = 0.0_wp
            return
        end if
        ! ln(lambda_D/b_min) for each field species, with Debye screening by
        ! electrons and both ions and b_min the larger of the classical and
        ! quantum impact parameters at vbar**2 = v0**2 + vth**2.
        debye = 1.0_wp/sqrt(4.0_wp*pi*e_cgs**2*(densi1*Z1**2/(tempi1*ev_cgs) &
                                 + densi2*Z2**2/(tempi2*ev_cgs) + dense/(tempe*ev_cgs)))
        field_mass = [am1*pmass, am2*pmass, emass]
        field_charge = [Z1, Z2, 1.0_wp]
        vth2 = [vti1**2, vti2**2, vte**2]
        do i = 1, nsorts
            reduced = mass*pmass*field_mass(i)/(mass*pmass + field_mass(i))
            vbar2 = v0**2 + vth2(i)
            bcl = charge*field_charge(i)*e_cgs**2/(reduced*vbar2)
            bqm = hbar/(2.0_wp*reduced*sqrt(vbar2))
            clog(i) = log(debye/max(bcl, bqm))
            bratio(i) = bcl/bqm
        end do
        if (any(clog <= 0.0_wp)) error stop &
            "Coulomb logarithm is not positive; plasma is not weakly coupled"
        alami1 = clog(1)
        alami2 = clog(2)
        alame = clog(3)
        if (present(clog_out)) clog_out = clog
        if (present(bratio_out)) bratio_out = bratio

        frecol_base = 2.d0*pi*dense*e_cgs**4*charge**2/((mass*pmass)**2*v0**3)
        frecol_base = frecol_base/v0

        efcolf_out(1) = frecol_base*Z1**2*alami1*densi1/dense
        efcolf_out(2) = frecol_base*Z2**2*alami2*densi2/dense
        efcolf_out(3) = frecol_base*alame

        efcolf_out = efcolf_out*velrat_out
    end subroutine compute_collision_coeffs

    subroutine compute_rates(v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        real(wp), intent(in) :: v0
        real(wp), intent(out) :: dchichi, slowrate, dchichi_norm, slowrate_norm

        real(wp) :: p, dpp, dhh, fpeff

        p = 1.d0
        call coleff(p, dpp, dhh, fpeff)

        dchichi = dhh*v0
        slowrate = abs(fpeff)*v0
        dchichi_norm = dhh
        slowrate_norm = abs(fpeff)
    end subroutine compute_rates

    subroutine init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0, &
                                       test_mass, test_charge)
        use simple_profiles, only: profile_max_s

        real(wp), intent(in) :: am1, am2, Z1, Z2, ealpha, v0
        real(wp), intent(in), optional :: test_mass, test_charge

        collision_s_max = profile_max_s()
        profile_am1 = am1
        profile_am2 = am2
        profile_Z1 = Z1
        profile_Z2 = Z2
        profile_ealpha = ealpha
        profile_test_mass = 4.0_wp
        profile_test_charge = 2.0_wp
        if (present(test_mass)) profile_test_mass = test_mass
        if (present(test_charge)) profile_test_charge = test_charge
        if (profile_test_mass <= 0.0_wp .or. profile_test_charge <= 0.0_wp) &
            error stop "Collision test species must have positive mass and charge"
        profile_ready = .true.
        uniform_ready = .false.
    end subroutine init_collision_profiles

    !> Without clog_in and bratio_in the Coulomb logarithm is speed independent.
    subroutine set_uniform_coeffs(efcolf_in, velrat_in, enrat_in, clog_in, bratio_in)
        real(wp), intent(in) :: efcolf_in(nsorts), velrat_in(nsorts), enrat_in(nsorts)
        real(wp), intent(in), optional :: clog_in(nsorts), bratio_in(nsorts)

        uniform_efcolf = efcolf_in
        uniform_velrat = velrat_in
        uniform_enrat = enrat_in
        uniform_clog = 0.0_wp
        uniform_bratio = 0.0_wp
        if (present(clog_in) .and. present(bratio_in)) then
            uniform_clog = clog_in
            uniform_bratio = bratio_in
        end if
        uniform_ready = .true.
    end subroutine set_uniform_coeffs

 subroutine get_local_coeffs(s, efcolf_loc, velrat_loc, enrat_loc, clog_loc, bratio_loc)
        use simple_profiles, only: get_plasma_params

        real(wp), intent(in) :: s
      real(wp), intent(out) :: efcolf_loc(nsorts), velrat_loc(nsorts), enrat_loc(nsorts)
        real(wp), intent(out), optional :: clog_loc(nsorts), bratio_loc(nsorts)

        real(wp) :: s_clamped, Te, Ti1, Ti2, ni1, ni2, v0_loc
        real(wp) :: densi1_loc, densi2_loc

        if (uniform_ready) then
            efcolf_loc = uniform_efcolf
            if (present(clog_loc)) clog_loc = uniform_clog
            if (present(bratio_loc)) bratio_loc = uniform_bratio
            velrat_loc = uniform_velrat
            enrat_loc = uniform_enrat
            return
        end if
        if (.not. profile_ready) then
            ! Without a plasma the operator is vacuum.
            efcolf_loc = 0.0_wp
            velrat_loc = 0.0_wp
            enrat_loc = 0.0_wp
            if (present(clog_loc)) clog_loc = 0.0_wp
            if (present(bratio_loc)) bratio_loc = 0.0_wp
            return
        end if
        ! Direct evaluation: a linear table in s misrepresents steep edge
        ! profiles by up to 89 % in pitch frequency with 100 intervals.
        s_clamped = max(0.0d0, min(collision_s_max, s))
        call get_plasma_params(s_clamped, Te, Ti1, Ti2, ni1, ni2)
        ! Preserve exact vacuum; retain the floor for positive densities.
        if (ni1 == 0.0_wp .and. ni2 == 0.0_wp) then
            densi1_loc = 0.0_wp
            densi2_loc = 0.0_wp
        else
            densi1_loc = max(ni1*1.0d-6, 1.0d0)
            densi2_loc = max(ni2*1.0d-6, 1.0d0)
        end if
        call compute_collision_coeffs(profile_am1, profile_am2, profile_Z1, &
                                      profile_Z2, densi1_loc, densi2_loc, &
                                      max(Ti1, 1.0d0), max(Ti2, 1.0d0), &
                                      max(Te, 1.0d0), profile_ealpha, v0_loc, &
                                      efcolf_loc, velrat_loc, enrat_loc, &
                                      profile_test_mass, profile_test_charge, &
                                      clog_loc, bratio_loc)
    end subroutine get_local_coeffs

    subroutine stost(z, dtauc, iswmode, ierr, max_substep_fraction, s_profile)
        use finite_check, only: is_finite

        integer, intent(in) :: iswmode
        real(wp), intent(in) :: dtauc
        real(wp), intent(inout) :: z(5)
        integer, intent(out) :: ierr
        real(wp), intent(in), optional :: max_substep_fraction
        !> Normalised toroidal flux of the profiles at z. Defaults to z(1),
        !> which is s only for integ_coords TEST, CANFLUX, VMEC and BOOZER;
        !> tracers pass profile_radius%profile_s(z(1:3)).
        real(wp), intent(in), optional :: s_profile

        real(wp), parameter :: default_fraction = 0.1_wp
        real(wp), parameter :: twopi = 8.0_wp*atan(1.0_wp)
        integer, parameter :: max_substeps = 100000
        real(wp) :: initial(5), fraction, remaining, step, p, dpp, dhh, fpeff
        real(wp) :: alam, angle, azimuth
        real(wp) :: p_switch, dpar, dperp, adrift, vhat(2), vnew(3), xi(3)
        real(wp) :: efcolf_loc(nsorts), velrat_loc(nsorts), enrat_loc(nsorts)
        real(wp) :: clog_loc(nsorts), bratio_loc(nsorts)
        real :: ur
        integer :: substep, i
        logical :: pitch, energy, diffusion, cartesian_ok

        initial = z
        ierr = 0
        fraction = default_fraction
        if (present(max_substep_fraction)) fraction = max_substep_fraction
        if (.not. all(is_finite(z)) .or. .not. is_finite(dtauc) .or. &
            .not. is_finite(fraction)) then
            ierr = 1
            return
        end if
        if (z(4) <= 0.0_wp .or. abs(z(5)) > 1.0_wp .or. dtauc < 0.0_wp .or. &
            fraction <= 0.0_wp .or. fraction > 0.25_wp .or. &
            iswmode < 1 .or. iswmode > 4) then
            ierr = 1
            return
        end if
        if (dtauc == 0.0_wp) return
        if (present(s_profile)) then
            if (.not. is_finite(s_profile)) then
                ierr = 1
                return
            end if
            call get_local_coeffs(s_profile, efcolf_loc, velrat_loc, enrat_loc, &
                                  clog_loc, bratio_loc)
        else
            call get_local_coeffs(z(1), efcolf_loc, velrat_loc, enrat_loc, &
                                  clog_loc, bratio_loc)
        end if
        ! A vacuum neither advances the state nor consumes random numbers.
        if (all(efcolf_loc == 0.0_wp)) return
        pitch = iswmode == 1 .or. iswmode == 4
        energy = iswmode /= 4
        diffusion = iswmode < 3
        remaining = dtauc
        ! Below the thermal speed of the slowest field species the radial
        ! process is Bessel-like and the Euler walk in p can collapse towards
        ! p = 0. There the substep is taken in Cartesian velocity space.
        cartesian_ok = iswmode <= 2 .and. any(efcolf_loc > 0.0_wp .and. &
                                              velrat_loc > 0.0_wp)
        p_switch = 0.0_wp
        if (cartesian_ok) p_switch = 1.0_wp/maxval(velrat_loc, &
                                                   mask=efcolf_loc > 0.0_wp)

        do substep = 1, max_substeps
            p = z(4)
            if (cartesian_ok .and. p < p_switch) then
                call coleff_cartesian(p, efcolf_loc, velrat_loc, enrat_loc, &
                                      dpar, dperp, adrift, clog_loc, bratio_loc)
                if (.not. all(is_finite([dpar, dperp, adrift])) .or. &
                    dpar < 0.0_wp .or. dperp < 0.0_wp) then
                    ierr = 1
                    exit
                end if
                step = remaining
                if (max(dpar, dperp) > 0.0_wp) step = min(step, &
                                      fraction**2*p_switch**2/(6.0_wp*max(dpar, dperp)))
                if (abs(adrift) > 0.0_wp) &
                    step = min(step, fraction*p_switch/abs(adrift))
                if (step <= 0.0_wp .or. remaining - step == remaining) then
                    ierr = 4 ! Cannot represent further time progress.
                    exit
                end if
                alam = z(5)
                vhat = [alam, sqrt(1.0_wp - alam**2)]
                do i = 1, 3
                    call getran(0, ur)
                    xi(i) = real(ur, wp)
                end do
                vnew(1:2) = p*vhat + (adrift*step + &
                                      sqrt(2.0_wp*dpar*step)*xi(1))*vhat + &
                            sqrt(2.0_wp*dperp*step)*xi(2)*[-vhat(2), vhat(1)]
                vnew(3) = sqrt(2.0_wp*dperp*step)*xi(3)
                z(4) = sqrt(sum(vnew**2))
                if (pitch .and. is_finite(z(4)) .and. z(4) > 0.0_wp) &
                    z(5) = max(-1.0_wp, min(1.0_wp, vnew(1)/z(4)))
                if (.not. all(is_finite(z)) .or. z(4) <= 0.0_wp .or. &
                    abs(z(5)) > 1.0_wp) then
                    ierr = 1
                    exit
                end if
                remaining = remaining - step
                if (remaining == 0.0_wp) return
                cycle
            end if
            call coleff_local(p, efcolf_loc, velrat_loc, enrat_loc, dpp, dhh, fpeff, &
                              clog_loc, bratio_loc)
            if (.not. all(is_finite([dpp, dhh, fpeff])) .or. &
                dpp < 0.0_wp .or. dhh < 0.0_wp) then
                ierr = 1
                exit
            end if
            ! Choose from the current state before drawing increments. Bound
            ! angular diffusion and both the relative drag and largest possible
            ! radial noise (getran(0) has support [-sqrt(3),sqrt(3)]).
            step = remaining
            if (pitch .and. dhh > 0.0_wp) &
                step = min(step, fraction**2/dhh)
            if (energy) then
                if (abs(fpeff) > 0.0_wp) &
                    step = min(step, fraction*p/abs(fpeff))
                if (diffusion .and. dpp > 0.0_wp) &
                    step = min(step, fraction**2*p**2/(6.0_wp*dpp))
            end if
            if (step <= 0.0_wp .or. remaining - step == remaining) then
                ierr = 4 ! Cannot represent further time progress.
                exit
            end if
            if (pitch) then
                ! Isotropic random rotation of fixed small angle. Its weak
                ! generator is dhh*d/dlambda[(1-lambda**2)*d/dlambda].
                ! Unlike Euler pitch kicks, it stays on the velocity sphere.
                angle = sqrt(4.0_wp*dhh*step)
                call random_number(ur)
                azimuth = twopi*real(ur, wp)
                alam = z(5)
                z(5) = alam*cos(angle) + sqrt(1.0_wp - alam**2)* &
                       sin(angle)*cos(azimuth)
            end if
            if (energy) then
                z(4) = p + fpeff*step
                if (diffusion) then
                    call getran(0, ur)
                    z(4) = z(4) + sqrt(2.0_wp*dpp*step)*real(ur, wp)
                end if
            end if
            if (.not. all(is_finite(z)) .or. z(4) <= 0.0_wp .or. &
                abs(z(5)) > 1.0_wp) then
                ierr = 1
                exit
            end if
            remaining = remaining - step
            if (remaining == 0.0_wp) return
        end do
        if (ierr == 0) ierr = 5 ! Explicit computational budget exhaustion.
        z = initial
    end subroutine stost

    subroutine getran(irand, ur)
        integer, intent(in) :: irand
        real, intent(out) :: ur

        call random_number(ur)

        if (irand == 0) then
            ur = 3.464102*(ur - 0.5)
        else
            if (ur > 0.5) then
                ur = 1.0
            else
                ur = -1.0
            end if
        end if
    end subroutine getran

end module collis_alp
