module magfie_sub
    use spline_vmec_sub, only: vmec_field
    use field_can_meiss, only: magfie_meiss
    use field_can_albert, only: magfie_albert
    use magfie_can_boozer_sub, only: magfie_can, magfie_boozer
    use util, only: twopi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use field_geoflux, only: geoflux_ready
    use geoflux_coordinates, only: geoflux_to_cyl
    use geoflux_field, only: splint_geoflux_field

    implicit none

! Define real(dp) kind parameter
    integer, parameter :: dp = kind(1.0d0)

    abstract interface
        subroutine magfie_base(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            import :: dp
            !            x(i)   - set of 3 curvilinear space coordinates (input)
            !            bmod   - dimensionless magnetic field module: bmod=B/B_ref
            !            sqrtg  - Jacobian of space coordinates (square root of
            !                     metric tensor
            !            bder   - derivatives of logarithm of bmod over space coords
            !                     (covariant vector)
            !            hcovar - covariant components of the unit vector along
            !                     the magnetic field
            !            hctrvr - contravariant components of the unit vector along
            !                     the magnetic field
            !            hcurl  - contravariant components of the curl of this vector
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: bmod, sqrtg
            real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
        end subroutine magfie_base
    end interface

    procedure(magfie_base), pointer :: magfie => null()

    integer, parameter :: TEST = -1, CANFLUX = 0, VMEC = 1, BOOZER = 2, MEISS = &
        3, ALBERT &
                          = 4, GEOFLUX = 5

contains

    subroutine init_magfie(id)
        integer, intent(in) :: id

        select case (id)
        case (TEST)
            magfie => magfie_test
        case (CANFLUX)
            magfie => magfie_can
        case (VMEC)
            if (geoflux_ready) then
                magfie => magfie_geoflux
            else
                magfie => magfie_vmec
            end if
        case (BOOZER)
            magfie => magfie_boozer
        case (MEISS)
            magfie => magfie_meiss
        case (ALBERT)
            magfie => magfie_albert
        case (GEOFLUX)
            magfie => magfie_geoflux
        case default
            print *, 'init_magfie: unknown id ', id
            error stop
        end select
    end subroutine init_magfie

    subroutine magfie_test(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        !> Magnetic field for analytic circular tokamak (TEST field).
        !>
        !> Coordinates: x(1)=s (flux-like), x(2)=theta (poloidal), x(3)=phi (toroidal).
        !> Mapping to minor radius: r = a*sqrt(s), with B0=1, R0=1, a=0.5, iota=1.
        !>
        !> This routine provides the full magfie interface (including hcurl) so that
        !> RK45 guiding-center integration has a consistent test field baseline.
        implicit none

        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(dp), parameter :: B0 = 1.0_dp, R0 = 1.0_dp, a = 0.5_dp, iota0 = 1.0_dp
        real(dp), parameter :: hs = 1.0d-4, ht = 1.0d-3*twopi
        real(dp) :: s, theta, phi, r, cth, sth, R_cyl
        real(dp) :: Ath, Aph, dAth_dr, dAph_dr
        real(dp) :: ds_fwd, ds_bwd, ds_den
        real(dp) :: bmod_plus, bmod_minus, bmod_tplus, bmod_tminus
        real(dp) :: hcov_plus(3), hcov_minus(3)
        real(dp) :: hcov_tplus(3), hcov_tminus(3)
        real(dp) :: dh_ds(3), dh_dt(3)
        real(dp) :: gss, gtt, gpp
        real(dp) :: Bsup_theta, Bsup_phi
        real(dp) :: Bcov_theta, Bcov_phi
        real(dp) :: sqrtg_geom

        s = max(0.0_dp, min(1.0_dp, x(1)))
        theta = x(2)
        phi = x(3)
        cth = cos(theta)
        sth = sin(theta)

        r = a*sqrt(s)
        R_cyl = R0 + r*cth

        ! Covariant vector potential (Ath, Aph) from field_can_test.
        Ath = B0*(r**2/2.0_dp - r**3/(3.0_dp*R0)*cth)
        Aph = -B0*iota0*(r**2/2.0_dp - r**4/(4.0_dp*a**2))

        ! Derivatives w.r.t r
        dAth_dr = B0*(r - r**2/R0*cth)
        dAph_dr = -B0*iota0*(r - r**3/a**2)

        ! Jacobian for (s,theta,phi) with r=a*sqrt(s):
        ! dV = (a^2/2) * R(s,theta) ds dtheta dphi.
        sqrtg_geom = 0.5_dp*a*a*R_cyl
        sqrtg = max(sqrtg_geom, 1.0d-14)

        ! Contravariant components of B = curl(A) in (s,theta,phi):
        ! B^s = (∂_θ A_φ - ∂_φ A_θ)/sqrtg = 0 for axisym + Aph(r).
        ! B^θ = -(∂_s A_φ)/sqrtg, B^φ = (∂_s A_θ)/sqrtg.
        if (s > 0.0_dp) then
            Bsup_theta = -(dAph_dr * (a/(2.0_dp*sqrt(s)))) / sqrtg
            Bsup_phi = (dAth_dr * (a/(2.0_dp*sqrt(s)))) / sqrtg
        else
            Bsup_theta = (B0*iota0)/max(R_cyl, 1.0d-12)
            Bsup_phi = B0/max(R_cyl, 1.0d-12)
        end if

        gss = (a*a)/(4.0_dp*max(s, 1.0d-14))
        gtt = r*r
        gpp = R_cyl*R_cyl

        Bcov_theta = gtt*Bsup_theta
        Bcov_phi = gpp*Bsup_phi

        bmod = sqrt((Bsup_theta*Bcov_theta) + (Bsup_phi*Bcov_phi))
        bmod = max(bmod, 1.0d-14)

        hcovar = 0.0_dp
        hcovar(2) = Bcov_theta/bmod
        hcovar(3) = Bcov_phi/bmod

        hctrvr = 0.0_dp
        hctrvr(2) = Bsup_theta/bmod
        hctrvr(3) = Bsup_phi/bmod

        ! Finite-difference derivatives for bder = ∂ ln(B)/∂x^i, and curl(h)^i.
        ds_fwd = min(hs, 1.0_dp - s)
        ds_bwd = min(hs, s)
        ds_den = ds_fwd + ds_bwd

        call magfie_test_eval_basic(s + ds_fwd, theta, phi, bmod_plus, hcov_plus)
        call magfie_test_eval_basic(s - ds_bwd, theta, phi, bmod_minus, hcov_minus)

        if (ds_den > 1.0d-16) then
            bder(1) = (bmod_plus - bmod_minus)/ds_den
            dh_ds = (hcov_plus - hcov_minus)/ds_den
        else
            bder(1) = 0.0_dp
            dh_ds = 0.0_dp
        end if

        call magfie_test_eval_basic(s, theta + ht, phi, bmod_tplus, hcov_tplus)
        call magfie_test_eval_basic(s, theta - ht, phi, bmod_tminus, hcov_tminus)
        bder(2) = (bmod_tplus - bmod_tminus)/(2.0_dp*ht)
        dh_dt = (hcov_tplus - hcov_tminus)/(2.0_dp*ht)

        bder(3) = 0.0_dp

        bder = bder/bmod

        ! curl(h)^s = (∂_θ h_φ - ∂_φ h_θ)/sqrtg, axisymmetric -> ∂_φ = 0.
        hcurl = 0.0_dp
        hcurl(1) = dh_dt(3)/sqrtg
        hcurl(2) = -dh_ds(3)/sqrtg
        hcurl(3) = dh_ds(2)/sqrtg

    end subroutine magfie_test

    subroutine magfie_test_eval_basic(s, theta, phi, bmod, hcov)
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: bmod, hcov(3)

        real(dp), parameter :: B0 = 1.0_dp, R0 = 1.0_dp, a = 0.5_dp, iota0 = 1.0_dp
        real(dp) :: s_clip, r, R_cyl, cth
        real(dp) :: Ath, Aph, dAth_dr, dAph_dr
        real(dp) :: sqrtg, gss, gtt, gpp
        real(dp) :: Bsup_theta, Bsup_phi, Bcov_theta, Bcov_phi

        s_clip = max(0.0_dp, min(1.0_dp, s))
        r = a*sqrt(s_clip)
        cth = cos(theta)
        R_cyl = R0 + r*cth

        Ath = B0*(r**2/2.0_dp - r**3/(3.0_dp*R0)*cth)
        Aph = -B0*iota0*(r**2/2.0_dp - r**4/(4.0_dp*a**2))
        dAth_dr = B0*(r - r**2/R0*cth)
        dAph_dr = -B0*iota0*(r - r**3/a**2)

        sqrtg = max(0.5_dp*a*a*R_cyl, 1.0d-14)

        if (s_clip > 0.0_dp) then
            Bsup_theta = -(dAph_dr * (a/(2.0_dp*sqrt(s_clip)))) / sqrtg
            Bsup_phi = (dAth_dr * (a/(2.0_dp*sqrt(s_clip)))) / sqrtg
        else
            Bsup_theta = (B0*iota0)/max(R_cyl, 1.0d-12)
            Bsup_phi = B0/max(R_cyl, 1.0d-12)
        end if

        gss = (a*a)/(4.0_dp*max(s_clip, 1.0d-14))
        gtt = r*r
        gpp = R_cyl*R_cyl

        Bcov_theta = gtt*Bsup_theta
        Bcov_phi = gpp*Bsup_phi

        bmod = sqrt((Bsup_theta*Bcov_theta) + (Bsup_phi*Bcov_phi))
        bmod = max(bmod, 1.0d-14)

        hcov = 0.0_dp
        hcov(2) = Bcov_theta/bmod
        hcov(3) = Bcov_phi/bmod
    end subroutine magfie_test_eval_basic

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine magfie_vmec(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        !
        ! Computes magnetic field module in units of the magnetic code  - bmod,
        ! square root of determinant of the metric tensor               - sqrtg,
        ! derivatives of the logarythm of the magnetic field module
        ! over coordinates                                              - bder,
        ! covariant componets of the unit vector of the magnetic
        ! field direction                                               - hcovar,
        ! contravariant components of this vector                       - hctrvr,
        ! contravariant component of the curl of this vector            - hcurl
        ! Order of coordinates is the following: x(1)=s (normalized toroidal flux),
        ! x(2)=theta (VMEC poloidal angle), x(3)=varphi (geometrical toroidal angle).
        !
        !  Input parameters:
        !            formal:  x(3)             - array of VMEC coordinates
        !  Output parameters:
        !            formal:  bmod
        !                     sqrtg
        !                     bder(3)          - derivatives of $\log(B)$
!                     hcovar(3)        - covariant components of unit vector $\bh$ along $\bB$
!                     hctrvr(3)        - contra-variant components of unit vector $\bh$ along $\bB$
     !                     hcurl(3)         - contra-variant components of curl of $\bh$
        !
        !  Called routines: vmec_field
        !
        implicit none
        !
        real(dp), parameter :: twopi = 2.d0*3.14159265358979d0, hs = 1.d-3, ht = &
            hs*twopi, &
                               hp = ht/5.d0
        !
        real(dp), intent(out) :: bmod, sqrtg
        real(dp) :: s, theta, varphi, A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                    sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                    Bcovar_r, Bcovar_vartheta, Bcovar_varphi
        real(dp) :: cjac, bcov_s_vmec, bcov_t_vmec, bcov_p_vmec
        real(dp) :: dhs_dt, dhs_dp, dht_ds, dht_dp, dhp_ds, dhp_dt
        real(dp), dimension(3), intent(in) :: x
        real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl
        !
        ! Begin derivatives over s
        !
        theta = x(2)
        varphi = x(3)
        s = x(1) + hs
        !
        call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, &
            dA_phi_ds, aiota, &
                        sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, &
                            Bctrvr_varphi, &
                        Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
        !
        bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
        bcov_s_vmec = Bcovar_r + dl_ds*Bcovar_vartheta
        bcov_t_vmec = (1.d0 + dl_dt)*Bcovar_vartheta
        bcov_p_vmec = Bcovar_varphi + dl_dp*Bcovar_vartheta
        bder(1) = bmod
        dht_ds = bcov_t_vmec/bmod
        dhp_ds = bcov_p_vmec/bmod
        !
        s = x(1) - hs
        !
        call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, &
            dA_phi_ds, aiota, &
                        sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, &
                            Bctrvr_varphi, &
                        Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
        !
        bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
        bcov_s_vmec = Bcovar_r + dl_ds*Bcovar_vartheta
        bcov_t_vmec = (1.d0 + dl_dt)*Bcovar_vartheta
        bcov_p_vmec = Bcovar_varphi + dl_dp*Bcovar_vartheta
        bder(1) = (bder(1) - bmod)/(2.d0*hs)
        dht_ds = (dht_ds - bcov_t_vmec/bmod)/(2.d0*hs)
        dhp_ds = (dhp_ds - bcov_p_vmec/bmod)/(2.d0*hs)
        !
        ! End derivatives over s
        !
        !-------------------------
        !
        ! Begin derivatives over theta
        !
        s = x(1)
        theta = x(2) + ht
        !
        call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, &
            dA_phi_ds, aiota, &
                        sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, &
                            Bctrvr_varphi, &
                        Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
        !
        bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
        bcov_s_vmec = Bcovar_r + dl_ds*Bcovar_vartheta
        bcov_t_vmec = (1.d0 + dl_dt)*Bcovar_vartheta
        bcov_p_vmec = Bcovar_varphi + dl_dp*Bcovar_vartheta
        bder(2) = bmod
        dhs_dt = bcov_s_vmec/bmod
        dhp_dt = bcov_p_vmec/bmod
        !
        theta = x(2) - ht
        !
        call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, &
            dA_phi_ds, aiota, &
                        sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, &
                            Bctrvr_varphi, &
                        Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
        !
        bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
        bcov_s_vmec = Bcovar_r + dl_ds*Bcovar_vartheta
        bcov_t_vmec = (1.d0 + dl_dt)*Bcovar_vartheta
        bcov_p_vmec = Bcovar_varphi + dl_dp*Bcovar_vartheta
        bder(2) = (bder(2) - bmod)/(2.d0*ht)
        dhs_dt = (dhs_dt - bcov_s_vmec/bmod)/(2.d0*ht)
        dhp_dt = (dhp_dt - bcov_p_vmec/bmod)/(2.d0*ht)
        !
        ! End derivatives over theta
        !
        !-------------------------
        !
        ! Begin derivatives over varphi
        !
        theta = x(2)
        varphi = x(3) + hp
        !
        call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, &
            dA_phi_ds, aiota, &
                        sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, &
                            Bctrvr_varphi, &
                        Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
        !
        bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
        bcov_s_vmec = Bcovar_r + dl_ds*Bcovar_vartheta
        bcov_t_vmec = (1.d0 + dl_dt)*Bcovar_vartheta
        bcov_p_vmec = Bcovar_varphi + dl_dp*Bcovar_vartheta
        bder(3) = bmod
        dhs_dp = bcov_s_vmec/bmod
        dht_dp = bcov_t_vmec/bmod
        !
        varphi = x(3) - hp
        !
        call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, &
            dA_phi_ds, aiota, &
                        sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, &
                            Bctrvr_varphi, &
                        Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
        !
        bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
        bcov_s_vmec = Bcovar_r + dl_ds*Bcovar_vartheta
        bcov_t_vmec = (1.d0 + dl_dt)*Bcovar_vartheta
        bcov_p_vmec = Bcovar_varphi + dl_dp*Bcovar_vartheta
        bder(3) = (bder(3) - bmod)/(2.d0*hp)
        dhs_dp = (dhs_dp - bcov_s_vmec/bmod)/(2.d0*hp)
        dht_dp = (dht_dp - bcov_t_vmec/bmod)/(2.d0*hp)
        !
        ! End derivatives over varphi
        !
        !-------------------------
        !
        varphi = x(3)
        !
        call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, &
            dA_phi_ds, aiota, &
                        sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, &
                            Bctrvr_varphi, &
                        Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
        !
        bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
        cjac = 1.d0 + dl_dt
        sqrtg = sqg*cjac
        bder = bder/bmod
        bcov_s_vmec = Bcovar_r + dl_ds*Bcovar_vartheta
        bcov_t_vmec = (1.d0 + dl_dt)*Bcovar_vartheta
        bcov_p_vmec = Bcovar_varphi + dl_dp*Bcovar_vartheta
        hcovar(1) = bcov_s_vmec/bmod
        hcovar(2) = bcov_t_vmec/bmod
        hcovar(3) = bcov_p_vmec/bmod
        hctrvr(1) = 0.d0
        hctrvr(2) = (Bctrvr_vartheta - dl_dp*Bctrvr_varphi)/(cjac*bmod)
        hctrvr(3) = Bctrvr_varphi/bmod
        hcurl(1) = (dhp_dt - dht_dp)/sqrtg
        hcurl(2) = (dhs_dp - dhp_ds)/sqrtg
        hcurl(3) = (dht_ds - dhs_dt)/sqrtg
        !
    end subroutine magfie_vmec
    !
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !

    subroutine magfie_geoflux(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(dp) :: s, theta, phi
        real(dp) :: ds_fwd, ds_bwd, ds_den
        real(dp) :: dt_step, dp_step
        real(dp) :: bmod_plus, bmod_minus
        real(dp) :: bmod_theta_plus, bmod_theta_minus
        real(dp) :: bmod_phi_plus, bmod_phi_minus
        real(dp) :: hcov_plus(3), hcov_minus(3)
        real(dp) :: hcov_theta_plus(3), hcov_theta_minus(3)
        real(dp) :: hcov_phi_plus(3), hcov_phi_minus(3)
        real(dp) :: basis(3, 3), g(3, 3), ginv(3, 3)
        real(dp) :: detg, sqrtg_geom
        real(dp) :: dh_ds(3), dh_dt(3), dh_dp(3)
        real(dp) :: phi_plus, phi_minus

        s = max(0.0_dp, min(1.0_dp, x(1)))
        theta = x(2)
        phi = x(3)

        call geoflux_eval_point(s, theta, phi, bmod, hcovar, sqrtg, basis, g, &
            ginv, detg, &
                                sqrtg_geom)

        if (sqrtg <= 0.0_dp) sqrtg = max(sqrtg_geom, 1.0d-12)
        sqrtg = max(sqrtg, 1.0d-12)

        if (.not. ieee_is_finite(bmod)) then
            error stop 'magfie_geoflux: non-finite Bmod'
        end if
        if (.not. all(ieee_is_finite(hcovar))) then
            error stop 'magfie_geoflux: non-finite hcovar'
        end if

        ds_fwd = min(1.0d-3, 1.0_dp - s)
        ds_bwd = min(1.0d-3, s)
        dt_step = 1.0d-3*twopi
        dp_step = dt_step/5.0d0

        call geoflux_eval_basic(s + ds_fwd, theta, phi, bmod_plus, hcov_plus)
        call geoflux_eval_basic(s - ds_bwd, theta, phi, bmod_minus, hcov_minus)

        ds_den = ds_fwd + ds_bwd
        if (ds_den > 1.0d-12) then
            bder(1) = (bmod_plus - bmod_minus)/ds_den
            dh_ds = (hcov_plus - hcov_minus)/ds_den
        else
            bder(1) = 0.0_dp
            dh_ds = 0.0_dp
        end if

        call geoflux_eval_basic(s, theta + dt_step, phi, bmod_theta_plus, &
            hcov_theta_plus)
        call geoflux_eval_basic(s, theta - dt_step, phi, bmod_theta_minus, &
            hcov_theta_minus)
        bder(2) = (bmod_theta_plus - bmod_theta_minus)/(2.0_dp*dt_step)
        dh_dt = (hcov_theta_plus - hcov_theta_minus)/(2.0_dp*dt_step)

        phi_plus = modulo(phi + dp_step, twopi)
        phi_minus = modulo(phi - dp_step, twopi)
        call geoflux_eval_basic(s, theta, phi_plus, bmod_phi_plus, hcov_phi_plus)
        call geoflux_eval_basic(s, theta, phi_minus, bmod_phi_minus, hcov_phi_minus)
        bder(3) = (bmod_phi_plus - bmod_phi_minus)/(2.0_dp*dp_step)
        dh_dp = (hcov_phi_plus - hcov_phi_minus)/(2.0_dp*dp_step)

        bder = bder/max(bmod, 1.0d-12)

    hctrvr = matmul(ginv, hcovar)

    if (sqrtg > 0.0_dp) then
      ! curl(h)^s = (∂_θ h_φ - ∂_φ h_θ) / sqrtg
      hcurl(1) = (dh_dt(3) - dh_dp(2))/sqrtg
      hcurl(2) = (dh_dp(1) - dh_ds(3))/sqrtg
      hcurl(3) = (dh_ds(2) - dh_dt(1))/sqrtg
    else
      hcurl = 0.0_dp
    end if

    end subroutine magfie_geoflux

    subroutine geoflux_eval_point(s, theta, phi, bmod, hcov, sqrtg, basis, g, ginv, &
                                  detg, sqrtg_geom)
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: bmod, hcov(3), sqrtg
        real(dp), intent(out) :: basis(3, 3), g(3, 3), ginv(3, 3)
        real(dp), intent(out) :: detg, sqrtg_geom
        real(dp) :: xcyl(3), jac(3, 3)
        real(dp) :: dRds, dZds, dRdtheta, dZdtheta, dRdphi, dZdphi
        real(dp) :: cosphi, sinphi
        real(dp) :: cross12(3)

        call geoflux_eval_basic(s, theta, phi, bmod, hcov, sqrtg, xcyl, jac)

        cosphi = cos(xcyl(2))
        sinphi = sin(xcyl(2))

        dRds = jac(1, 1)
        dZds = jac(3, 1)
        dRdtheta = jac(1, 2)
        dZdtheta = jac(3, 2)
        dRdphi = jac(1, 3)
        dZdphi = jac(3, 3)

        basis(:, 1) = (/dRds*cosphi, dRds*sinphi, dZds/)
        basis(:, 2) = (/dRdtheta*cosphi, dRdtheta*sinphi, dZdtheta/)
        basis(:, 3) = (/dRdphi*cosphi - xcyl(1)*sinphi, &
                        dRdphi*sinphi + xcyl(1)*cosphi, dZdphi/)

        call compute_metric(basis, g, ginv, detg)
        call cross_product(basis(:, 2), basis(:, 3), cross12)
        sqrtg_geom = abs(dot_product(basis(:, 1), cross12))
        sqrtg = max(sqrtg, sqrtg_geom)
    end subroutine geoflux_eval_point

    subroutine geoflux_eval_basic(s, theta, phi, bmod, hcov, sqrtg, xcyl, jac)
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: bmod, hcov(3)
        real(dp), intent(out), optional :: sqrtg
        real(dp), intent(out), optional :: xcyl(3), jac(3, 3)

        real(dp) :: s_clip
        real(dp) :: sqg_tmp(3)
        real(dp) :: acov_tmp(3), hcov_tmp(3)
        real(dp) :: cyl_tmp(3), jac_tmp(3, 3)

        s_clip = max(0.0_dp, min(1.0_dp, s))

        call geoflux_to_cyl((/s_clip, theta, phi/), cyl_tmp, jac_tmp)
        call splint_geoflux_field(s_clip, theta, phi, acov_tmp, hcov_tmp, bmod, sqg_tmp)

        hcov = hcov_tmp

        if (present(sqrtg)) sqrtg = abs(sqg_tmp(1))
        if (present(xcyl)) xcyl = cyl_tmp
        if (present(jac)) jac = jac_tmp
    end subroutine geoflux_eval_basic

    subroutine compute_metric(basis, g, ginv, detg)
        real(dp), intent(in) :: basis(3, 3)
        real(dp), intent(out) :: g(3, 3), ginv(3, 3)
        real(dp), intent(out) :: detg
        integer :: i, j

        do i = 1, 3
            do j = 1, 3
                g(i, j) = dot_product(basis(:, i), basis(:, j))
            end do
        end do

        call invert3x3(g, ginv, detg)
        if (abs(detg) < 1.0d-16) then
            detg = 1.0d0
            ginv = 0.0_dp
            do i = 1, 3
                ginv(i, i) = 1.0_dp
            end do
        end if
    end subroutine compute_metric

    subroutine cross_product(a, b, c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp), intent(out) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end subroutine cross_product

    subroutine invert3x3(a, ainv, det)
        real(dp), intent(in) :: a(3, 3)
        real(dp), intent(out) :: ainv(3, 3)
        real(dp), intent(out) :: det

        det = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) &
              - a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) &
              + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))

        if (abs(det) < 1.0d-16) then
            det = 0.0_dp
            ainv = 0.0_dp
            return
        end if

        ainv(1, 1) = (a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2))/det
        ainv(1, 2) = -(a(1, 2)*a(3, 3) - a(1, 3)*a(3, 2))/det
        ainv(1, 3) = (a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2))/det
        ainv(2, 1) = -(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1))/det
        ainv(2, 2) = (a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1))/det
        ainv(2, 3) = -(a(1, 1)*a(2, 3) - a(1, 3)*a(2, 1))/det
        ainv(3, 1) = (a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))/det
        ainv(3, 2) = -(a(1, 1)*a(3, 2) - a(1, 2)*a(3, 1))/det
        ainv(3, 3) = (a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1))/det
    end subroutine invert3x3

end module magfie_sub
