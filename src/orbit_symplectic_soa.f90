module orbit_symplectic_soa
use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_boozer, only: eval_field_booz_many
#ifdef SIMPLE_PROFILE_COUNTERS
use profile_counters, only: prof_add_soa_newton_call, prof_add_soa_newton_iter, &
    prof_add_soa_eval_newton
#endif
use vector_potentail_mod, only: torflux

implicit none

#ifdef BATCH_SIZE
integer, parameter :: BATCH_PTS = BATCH_SIZE
#else
integer, parameter :: BATCH_PTS = 256
#endif

contains

subroutine get_derivatives_many(npts, ro0, mu, pphi, &
        Ath, Aph, dAth_dr, dAph_dr, hth, hph, dhth, dhph, Bmod, dBmod, &
        vpar, dvpar, pth, dpth, H, dH)
    integer, intent(in) :: npts
    real(dp), intent(in) :: ro0
    real(dp), intent(in) :: mu(npts)
    real(dp), intent(in) :: pphi(npts)
    real(dp), intent(in) :: Ath(npts), Aph(npts)
    real(dp), intent(in) :: dAth_dr(npts), dAph_dr(npts)
    real(dp), intent(in) :: hth(npts), hph(npts)
    real(dp), intent(in) :: dhth(3, npts), dhph(3, npts)
    real(dp), intent(in) :: Bmod(npts), dBmod(3, npts)
    real(dp), intent(out) :: vpar(npts), dvpar(4, npts)
    real(dp), intent(out) :: pth(npts), dpth(4, npts)
    real(dp), intent(out) :: H(npts), dH(4, npts)

    integer :: i

    !$acc kernels present(mu, pphi, Ath, Aph, dAth_dr, dAph_dr, hth, hph, &
    !$acc&                dhth, dhph, Bmod, dBmod, vpar, dvpar, pth, dpth, H, dH)
    !$omp simd
    do i = 1, npts
        vpar(i) = (pphi(i) - Aph(i) / ro0) / hph(i)
        dvpar(1, i) = (-dAph_dr(i) / ro0 - dhph(1, i) * vpar(i)) / hph(i)
        dvpar(2, i) = -dhph(2, i) * vpar(i) / hph(i)
        dvpar(3, i) = -dhph(3, i) * vpar(i) / hph(i)
        dvpar(4, i) = 1.0d0 / hph(i)

        pth(i) = vpar(i) * hth(i) + Ath(i) / ro0
        dpth(1, i) = dvpar(1, i) * hth(i) + vpar(i) * dhth(1, i) + dAth_dr(i) / ro0
        dpth(2, i) = dvpar(2, i) * hth(i) + vpar(i) * dhth(2, i)
        dpth(3, i) = dvpar(3, i) * hth(i) + vpar(i) * dhth(3, i)
        dpth(4, i) = dvpar(4, i) * hth(i)

        H(i) = 0.5d0 * vpar(i)**2 + mu(i) * Bmod(i)
        dH(1, i) = vpar(i) * dvpar(1, i) + mu(i) * dBmod(1, i)
        dH(2, i) = vpar(i) * dvpar(2, i) + mu(i) * dBmod(2, i)
        dH(3, i) = vpar(i) * dvpar(3, i) + mu(i) * dBmod(3, i)
        dH(4, i) = vpar(i) * dvpar(4, i)
    end do
    !$omp end simd
    !$acc end kernels
end subroutine get_derivatives_many


subroutine get_derivatives2_many(npts, ro0, mu, pphi, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod, &
        vpar, dvpar, d2vpar, pth, dpth, d2pth, H, dH, d2H)
    integer, intent(in) :: npts
    real(dp), intent(in) :: ro0
    real(dp), intent(in) :: mu(npts)
    real(dp), intent(in) :: pphi(npts)
    real(dp), intent(in) :: Ath(npts), Aph(npts)
    real(dp), intent(in) :: dAth_dr(npts), dAph_dr(npts), d2Aph_dr2(npts)
    real(dp), intent(in) :: hth(npts), hph(npts)
    real(dp), intent(in) :: dhth(3, npts), dhph(3, npts)
    real(dp), intent(in) :: d2hth_dr2(npts), d2hph_dr2(npts)
    real(dp), intent(in) :: Bmod(npts), dBmod(3, npts), d2Bmod(6, npts)
    real(dp), intent(out) :: vpar(npts), dvpar(4, npts), d2vpar(10, npts)
    real(dp), intent(out) :: pth(npts), dpth(4, npts), d2pth(10, npts)
    real(dp), intent(out) :: H(npts), dH(4, npts), d2H(10, npts)

    integer :: i

    call get_derivatives_many(npts, ro0, mu, pphi, &
        Ath, Aph, dAth_dr, dAph_dr, hth, hph, dhth, dhph, Bmod, dBmod, &
        vpar, dvpar, pth, dpth, H, dH)

    !$acc kernels
    !$omp simd
    do i = 1, npts
        d2vpar(1, i) = (-d2Aph_dr2(i) / ro0 - d2hph_dr2(i) * vpar(i)) / hph(i) &
            - 2.0d0 * dhph(1, i) * dvpar(1, i) / hph(i)
        d2vpar(7, i) = -dhph(1, i) * dvpar(4, i) / hph(i)

        d2H(1, i) = vpar(i) * d2vpar(1, i) + dvpar(1, i)**2 + mu(i) * d2Bmod(1, i)
        d2H(2, i) = vpar(i) * dvpar(1, i) * dvpar(2, i) / vpar(i) + mu(i) * d2Bmod(2, i)
        d2H(3, i) = vpar(i) * dvpar(1, i) * dvpar(3, i) / vpar(i) + mu(i) * d2Bmod(3, i)
        d2H(7, i) = vpar(i) * d2vpar(7, i) + dvpar(1, i) * dvpar(4, i)
        d2H(8, i) = dvpar(2, i) * dvpar(4, i)
        d2H(9, i) = dvpar(3, i) * dvpar(4, i)

        d2pth(1, i) = d2vpar(1, i) * hth(i) + 2.0d0 * dvpar(1, i) * dhth(1, i) &
            + vpar(i) * d2hth_dr2(i)
        d2pth(2, i) = dvpar(1, i) * dhth(2, i) + dvpar(2, i) * dhth(1, i)
        d2pth(3, i) = dvpar(1, i) * dhth(3, i) + dvpar(3, i) * dhth(1, i)
        d2pth(7, i) = d2vpar(7, i) * hth(i) + dvpar(4, i) * dhth(1, i)
        d2pth(8, i) = dvpar(4, i) * dhth(2, i)
        d2pth(9, i) = dvpar(4, i) * dhth(3, i)
    end do
    !$omp end simd
    !$acc end kernels
end subroutine get_derivatives2_many


subroutine jac_sympl_euler1_many(npts, dt, z_pphi, pthold, x_pphi, fjac, &
        pth, dpth, d2pth, dH, d2H)
    !> Compute Jacobian using precomputed values (get_derivatives2_many).
    integer, intent(in) :: npts
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: z_pphi(npts), pthold(npts)
    real(dp), intent(in) :: x_pphi(npts)
    real(dp), intent(out) :: fjac(2, 2, npts)
    real(dp), intent(in) :: pth(npts), dpth(4, npts), d2pth(10, npts)
    real(dp), intent(in) :: dH(4, npts), d2H(10, npts)

    integer :: i

    !$acc kernels
    !$omp simd
    do i = 1, npts
        fjac(1, 1, i) = d2pth(1, i) * (pth(i) - pthold(i)) + dpth(1, i)**2 &
            + dt * (d2H(2, i) * dpth(1, i) + dH(2, i) * d2pth(1, i) &
                  - d2H(1, i) * dpth(2, i) - dH(1, i) * d2pth(2, i))
        fjac(1, 2, i) = d2pth(7, i) * (pth(i) - pthold(i)) + dpth(1, i) * dpth(4, i) &
            + dt * (d2H(8, i) * dpth(1, i) + dH(2, i) * d2pth(7, i) &
                  - d2H(7, i) * dpth(2, i) - dH(1, i) * d2pth(8, i))
        fjac(2, 1, i) = d2pth(1, i) * (x_pphi(i) - z_pphi(i)) &
            + dt * (d2H(3, i) * dpth(1, i) + dH(3, i) * d2pth(1, i) &
                  - d2H(1, i) * dpth(3, i) - dH(1, i) * d2pth(3, i))
        fjac(2, 2, i) = d2pth(7, i) * (x_pphi(i) - z_pphi(i)) + dpth(1, i) &
            + dt * (d2H(9, i) * dpth(1, i) + dH(3, i) * d2pth(7, i) &
                  - d2H(7, i) * dpth(3, i) - dH(1, i) * d2pth(9, i))
    end do
    !$omp end simd
    !$acc end kernels
end subroutine jac_sympl_euler1_many


subroutine newton1_soa(npts, dt, ro0, mu, atol, rtol, maxit, &
        z_r, z_th, z_ph, z_pphi, pthold_in, x_r, x_pphi, converged, &
        xlast_r_out, xlast_pphi_out, &
        pth_out, dpth_out, dH_out, d2pth_1_out, d2pth_7_out, d2H_1_out, d2H_7_out, &
        vpar_out, dvpar_out, hth_out, hph_out, dhth_out, dhph_out)
    integer, intent(in) :: npts, maxit
    real(dp), intent(in) :: dt, ro0, atol, rtol
    real(dp), intent(in) :: mu(npts)
    real(dp), intent(in) :: z_r(npts), z_th(npts), z_ph(npts), z_pphi(npts)
    real(dp), intent(in), optional :: pthold_in(npts)
    real(dp), intent(inout) :: x_r(npts), x_pphi(npts)
    logical, intent(out) :: converged(npts)
    real(dp), intent(out) :: xlast_r_out(npts), xlast_pphi_out(npts)
    real(dp), intent(out) :: pth_out(npts), dpth_out(4, npts)
    real(dp), intent(out) :: dH_out(4, npts)
    real(dp), intent(out) :: d2pth_1_out(npts), d2pth_7_out(npts)
    real(dp), intent(out) :: d2H_1_out(npts), d2H_7_out(npts)
    real(dp), intent(out) :: vpar_out(npts), dvpar_out(4, npts)
    real(dp), intent(out) :: hth_out(npts), hph_out(npts)
    real(dp), intent(out) :: dhth_out(3, npts), dhph_out(3, npts)

    real(dp) :: Ath(BATCH_PTS), Aph(BATCH_PTS), dAth_dr(BATCH_PTS)
    real(dp) :: dAph_dr(BATCH_PTS), d2Aph_dr2(BATCH_PTS)
    real(dp) :: hth(BATCH_PTS), hph(BATCH_PTS)
    real(dp) :: dhth(3, BATCH_PTS), dhph(3, BATCH_PTS)
    real(dp) :: d2hth_dr2(BATCH_PTS), d2hph_dr2(BATCH_PTS)
    real(dp) :: Bmod(BATCH_PTS), dBmod(3, BATCH_PTS), d2Bmod(6, BATCH_PTS)
    real(dp) :: pth(BATCH_PTS), dpth(4, BATCH_PTS), H(BATCH_PTS), dH(4, BATCH_PTS)
    real(dp) :: d2pth(10, BATCH_PTS), d2H(10, BATCH_PTS)
    real(dp) :: vpar(BATCH_PTS), dvpar(4, BATCH_PTS)
    real(dp) :: fvec(2, BATCH_PTS), fjac(2, 2, BATCH_PTS)
    real(dp) :: xlast_r(BATCH_PTS), xlast_pphi(BATCH_PTS)
    real(dp) :: pthold_local(BATCH_PTS)
    real(dp) :: d2vpar(10, BATCH_PTS)
    real(dp) :: tolref_pphi(BATCH_PTS)
    logical :: just_converged
    real(dp) :: det, ijac11, ijac12, ijac21, ijac22
    real(dp) :: dx_r, dx_pphi
    integer :: kit, i, j

    converged(1:npts) = .false.
    tolref_pphi(1:npts) = abs(10.0d0 * torflux / ro0)
    if (present(pthold_in)) pthold_local(1:npts) = pthold_in(1:npts)

#ifdef SIMPLE_PROFILE_COUNTERS
    call prof_add_soa_newton_call
#endif
    do kit = 1, maxit
#ifdef SIMPLE_PROFILE_COUNTERS
        call prof_add_soa_newton_iter(kit)
        call prof_add_soa_eval_newton(npts)
#endif
        call eval_field_booz_many(npts, x_r, z_th, z_ph, &
            Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
            hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
            Bmod, dBmod, d2Bmod)

        call get_derivatives2_many(npts, ro0, mu, x_pphi, &
            Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
            hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
            Bmod, dBmod, d2Bmod, &
            vpar, dvpar, d2vpar, pth, dpth, d2pth, H, dH, d2H)

        if ((.not. present(pthold_in)) .and. kit == 1) then
            pthold_local(1:npts) = pth(1:npts)
        end if

        !$acc kernels
        !$omp simd
        do i = 1, npts
            fvec(1, i) = dpth(1, i) * (pth(i) - pthold_local(i)) &
                + dt * (dH(2, i) * dpth(1, i) - dH(1, i) * dpth(2, i))
            fvec(2, i) = dpth(1, i) * (x_pphi(i) - z_pphi(i)) &
                + dt * (dH(3, i) * dpth(1, i) - dH(1, i) * dpth(3, i))
        end do
        !$omp end simd
        !$acc end kernels

        call jac_sympl_euler1_many(npts, dt, z_pphi, pthold_local, x_pphi, fjac, &
            pth, dpth, d2pth, dH, d2H)

        do i = 1, npts
            if (converged(i)) cycle

            xlast_r(i) = x_r(i)
            xlast_pphi(i) = x_pphi(i)

            det = fjac(1, 1, i) * fjac(2, 2, i) - fjac(1, 2, i) * fjac(2, 1, i)
            ijac11 = fjac(2, 2, i) / det
            ijac12 = -fjac(1, 2, i) / det
            ijac21 = -fjac(2, 1, i) / det
            ijac22 = fjac(1, 1, i) / det

            dx_r = ijac11 * fvec(1, i) + ijac12 * fvec(2, i)
            dx_pphi = ijac21 * fvec(1, i) + ijac22 * fvec(2, i)

            x_r(i) = x_r(i) - dx_r
            x_pphi(i) = x_pphi(i) - dx_pphi

            if (x_r(i) < 0.0d0) x_r(i) = 0.01d0

            tolref_pphi(i) = max(abs(x_pphi(i)), tolref_pphi(i))

            just_converged = .false.
            if (abs(fvec(1, i)) < atol .and. abs(fvec(2, i)) < atol) then
                just_converged = .true.
            else if (abs(x_r(i) - xlast_r(i)) < rtol .and. &
                     abs(x_pphi(i) - xlast_pphi(i)) < rtol * tolref_pphi(i)) then
                just_converged = .true.
            end if

            if (just_converged) then
                converged(i) = .true.
                xlast_r_out(i) = xlast_r(i)
                xlast_pphi_out(i) = xlast_pphi(i)
                pth_out(i) = pth(i)
                do j = 1, 4
                    dpth_out(j, i) = dpth(j, i)
                    dH_out(j, i) = dH(j, i)
                    dvpar_out(j, i) = dvpar(j, i)
                end do
                d2pth_1_out(i) = d2pth(1, i)
                d2pth_7_out(i) = d2pth(7, i)
                d2H_1_out(i) = d2H(1, i)
                d2H_7_out(i) = d2H(7, i)
                vpar_out(i) = vpar(i)
                hth_out(i) = hth(i)
                hph_out(i) = hph(i)
                do j = 1, 3
                    dhth_out(j, i) = dhth(j, i)
                    dhph_out(j, i) = dhph(j, i)
                end do
            end if
        end do

        if (all(converged)) exit
    end do

    do i = 1, npts
        if (.not. converged(i)) then
            xlast_r_out(i) = xlast_r(i)
            xlast_pphi_out(i) = xlast_pphi(i)
            pth_out(i) = pth(i)
            do j = 1, 4
                dpth_out(j, i) = dpth(j, i)
                dH_out(j, i) = dH(j, i)
                dvpar_out(j, i) = dvpar(j, i)
            end do
            d2pth_1_out(i) = d2pth(1, i)
            d2pth_7_out(i) = d2pth(7, i)
            d2H_1_out(i) = d2H(1, i)
            d2H_7_out(i) = d2H(7, i)
            vpar_out(i) = vpar(i)
            hth_out(i) = hth(i)
            hph_out(i) = hph(i)
            do j = 1, 3
                dhth_out(j, i) = dhth(j, i)
                dhph_out(j, i) = dhph(j, i)
            end do
        end if
    end do
end subroutine newton1_soa


subroutine orbit_timestep_euler1_soa(npts, dt, ntau, ro0, mu, atol, rtol, maxit, &
        z_r, z_th, z_ph, z_pphi, escaped, ierr)
    integer, intent(in) :: npts, ntau, maxit
    real(dp), intent(in) :: dt, ro0, atol, rtol
    real(dp), intent(in) :: mu(npts)
    real(dp), intent(inout) :: z_r(npts), z_th(npts), z_ph(npts), z_pphi(npts)
    logical, intent(out) :: escaped(npts)
    integer, intent(out) :: ierr(npts)

    real(dp) :: pthold(BATCH_PTS), x_r(BATCH_PTS), x_pphi(BATCH_PTS)
    real(dp) :: xlast_r(BATCH_PTS), xlast_pphi(BATCH_PTS)
    real(dp) :: hth(BATCH_PTS), hph(BATCH_PTS)
    real(dp) :: dhth(3, BATCH_PTS), dhph(3, BATCH_PTS)
    real(dp) :: vpar(BATCH_PTS), dvpar(4, BATCH_PTS)
    real(dp) :: pth(BATCH_PTS), dpth(4, BATCH_PTS), dH(4, BATCH_PTS)
    real(dp) :: d2pth_1(BATCH_PTS), d2pth_7(BATCH_PTS)
    real(dp) :: d2H_1(BATCH_PTS), d2H_7(BATCH_PTS)
    logical :: converged(BATCH_PTS)
    real(dp) :: dr, dpphi
    integer :: ktau, i

    escaped = .false.
    ierr = 0

    do ktau = 1, ntau
        x_r(1:npts) = z_r(1:npts)
        x_pphi(1:npts) = z_pphi(1:npts)

        if (ktau > 1) pthold(1:npts) = pth(1:npts)

        if (ktau == 1) then
            call newton1_soa(npts, dt, ro0, mu, atol, rtol, maxit, &
                z_r, z_th, z_ph, z_pphi, x_r=x_r, x_pphi=x_pphi, converged=converged, &
                xlast_r_out=xlast_r, xlast_pphi_out=xlast_pphi, &
                pth_out=pth, dpth_out=dpth, dH_out=dH, &
                d2pth_1_out=d2pth_1, d2pth_7_out=d2pth_7, d2H_1_out=d2H_1, d2H_7_out=d2H_7, &
                vpar_out=vpar, dvpar_out=dvpar, hth_out=hth, hph_out=hph, dhth_out=dhth, dhph_out=dhph)
        else
            call newton1_soa(npts, dt, ro0, mu, atol, rtol, maxit, &
                z_r, z_th, z_ph, z_pphi, pthold_in=pthold, x_r=x_r, x_pphi=x_pphi, converged=converged, &
                xlast_r_out=xlast_r, xlast_pphi_out=xlast_pphi, &
                pth_out=pth, dpth_out=dpth, dH_out=dH, &
                d2pth_1_out=d2pth_1, d2pth_7_out=d2pth_7, d2H_1_out=d2H_1, d2H_7_out=d2H_7, &
                vpar_out=vpar, dvpar_out=dvpar, hth_out=hth, hph_out=hph, dhth_out=dhth, dhph_out=dhph)
        end if

        !$omp simd private(dr, dpphi)
        do i = 1, npts
            if (escaped(i)) cycle

            if (x_r(i) > 1.0d0) then
                escaped(i) = .true.
                ierr(i) = 1
                cycle
            end if

            if (x_r(i) < 0.0d0) x_r(i) = 0.01d0

            z_r(i) = x_r(i)
            z_pphi(i) = x_pphi(i)

            dr = x_r(i) - xlast_r(i)
            dpphi = x_pphi(i) - xlast_pphi(i)
            pth(i) = pth(i) + dpth(1, i) * dr + dpth(4, i) * dpphi
            dH(1, i) = dH(1, i) + d2H_1(i) * dr + d2H_7(i) * dpphi
            dpth(1, i) = dpth(1, i) + d2pth_1(i) * dr + d2pth_7(i) * dpphi
            vpar(i) = vpar(i) + dvpar(1, i) * dr + dvpar(4, i) * dpphi
            hth(i) = hth(i) + dhth(1, i) * dr
            hph(i) = hph(i) + dhph(1, i) * dr
        end do
        !$omp end simd

        !$omp simd
        do i = 1, npts
            if (escaped(i)) cycle
            z_th(i) = z_th(i) + dt * dH(1, i) / dpth(1, i)
            z_ph(i) = z_ph(i) + dt * (vpar(i) - dH(1, i) / dpth(1, i) * hth(i)) / hph(i)
        end do
        !$omp end simd
    end do
end subroutine orbit_timestep_euler1_soa


subroutine trace_orbit_soa(npts, zstart, ntimstep, ntau, dt_norm, ro0_in, &
        atol, rtol, maxit, z_final, times_lost, ierr)
    use boozer_sub, only: vmec_to_boozer
    integer, intent(in) :: npts, ntimstep, ntau, maxit
    real(dp), intent(in) :: zstart(5, npts)
    real(dp), intent(in) :: dt_norm, ro0_in, atol, rtol
    real(dp), intent(out) :: z_final(5, npts)
    real(dp), intent(out) :: times_lost(npts)
    integer, intent(out) :: ierr(npts)

    real(dp) :: z_r(BATCH_PTS), z_th(BATCH_PTS), z_ph(BATCH_PTS), z_pphi(BATCH_PTS)
    real(dp) :: mu(BATCH_PTS), pabs(BATCH_PTS), vpar(BATCH_PTS)
    real(dp) :: Ath(BATCH_PTS), Aph(BATCH_PTS), dAth_dr(BATCH_PTS)
    real(dp) :: dAph_dr(BATCH_PTS), d2Aph_dr2(BATCH_PTS)
    real(dp) :: hth(BATCH_PTS), hph(BATCH_PTS)
    real(dp) :: dhth(3, BATCH_PTS), dhph(3, BATCH_PTS)
    real(dp) :: d2hth_dr2(BATCH_PTS), d2hph_dr2(BATCH_PTS)
    real(dp) :: Bmod(BATCH_PTS), dBmod(3, BATCH_PTS), d2Bmod(6, BATCH_PTS)
    logical :: escaped(BATCH_PTS)
    real(dp) :: ro0_norm, dt, theta_B, phi_B
    integer :: i, it

    ro0_norm = ro0_in / dsqrt(2.0d0)
    dt = dt_norm / dsqrt(2.0d0)

    do i = 1, npts
        z_r(i) = zstart(1, i)
        call vmec_to_boozer(zstart(1, i), mod(zstart(2, i), 6.283185307179586d0), &
            mod(zstart(3, i), 6.283185307179586d0), theta_B, phi_B)
        z_th(i) = mod(theta_B, 6.283185307179586d0)
        z_ph(i) = mod(phi_B, 6.283185307179586d0)
        pabs(i) = zstart(4, i)
    end do

    call eval_field_booz_many(npts, z_r, z_th, z_ph, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod)

    !$acc kernels
    !$omp simd
    do i = 1, npts
        mu(i) = 0.5d0 * pabs(i)**2 * (1.0d0 - zstart(5, i)**2) / Bmod(i) * 2.0d0
        vpar(i) = pabs(i) * zstart(5, i) * dsqrt(2.0d0)
        z_pphi(i) = vpar(i) * hph(i) + Aph(i) / ro0_norm
    end do
    !$omp end simd
    !$acc end kernels

    escaped(1:npts) = .false.
    ierr(1:npts) = 0
    times_lost(1:npts) = 0.0d0

    do it = 1, ntimstep
        call orbit_timestep_euler1_soa(npts, dt, ntau, ro0_norm, mu, &
            atol, rtol, maxit, z_r, z_th, z_ph, z_pphi, escaped, ierr)

        !$omp simd
        do i = 1, npts
            if (escaped(i) .and. times_lost(i) < 1.0d-30) then
                times_lost(i) = real(it * ntau, dp) * dt_norm
            end if
        end do
        !$omp end simd

        if (all(escaped(1:npts))) exit
    end do

    call eval_field_booz_many(npts, z_r, z_th, z_ph, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod)

    !$acc kernels
    !$omp simd
    do i = 1, npts
        z_final(1, i) = z_r(i)
        z_final(2, i) = z_th(i)
        z_final(3, i) = z_ph(i)
        vpar(i) = (z_pphi(i) - Aph(i) / ro0_norm) / hph(i)
        z_final(4, i) = dsqrt(mu(i) * Bmod(i) + 0.5d0 * vpar(i)**2)
        z_final(5, i) = vpar(i) / (z_final(4, i) * dsqrt(2.0d0))
    end do
    !$omp end simd
    !$acc end kernels
end subroutine trace_orbit_soa


subroutine trace_orbit_soa_omp(npts, zstart, ntimstep, ntau, dt_norm, ro0_in, &
        atol, rtol, maxit, z_final, times_lost, ierr)
    use boozer_sub, only: vmec_to_boozer
    integer, intent(in) :: npts, ntimstep, ntau, maxit
    real(dp), intent(in) :: zstart(5, npts)
    real(dp), intent(in) :: dt_norm, ro0_in, atol, rtol
    real(dp), intent(out) :: z_final(5, npts)
    real(dp), intent(out) :: times_lost(npts)
    integer, intent(out) :: ierr(npts)

    integer :: istart, iend, npts_chunk, thread_id, nthreads
    integer :: omp_get_thread_num, omp_get_num_threads

    !$omp parallel private(istart, iend, npts_chunk, thread_id) &
    !$omp& shared(nthreads)
    thread_id = 0
    nthreads = 1
    !$ thread_id = omp_get_thread_num()
    !$ nthreads = omp_get_num_threads()

    npts_chunk = (npts + nthreads - 1) / nthreads
    istart = thread_id * npts_chunk + 1
    iend = min(istart + npts_chunk - 1, npts)

    if (istart <= npts) then
        call trace_orbit_soa_chunk(iend - istart + 1, &
            zstart(:, istart:iend), ntimstep, ntau, dt_norm, ro0_in, &
            atol, rtol, maxit, z_final(:, istart:iend), &
            times_lost(istart:iend), ierr(istart:iend))
    end if
    !$omp end parallel
end subroutine trace_orbit_soa_omp


subroutine trace_orbit_soa_omp1(npts, zstart, ntimstep, ntau, dt_norm, ro0_in, &
        atol, rtol, maxit, z_final, times_lost, ierr)
    !> OMP parallel do over particles with batch size 1.
    !> Each thread processes exactly 1 particle - equivalent memory access to AoS.
    integer, intent(in) :: npts, ntimstep, ntau, maxit
    real(dp), intent(in) :: zstart(5, npts)
    real(dp), intent(in) :: dt_norm, ro0_in, atol, rtol
    real(dp), intent(out) :: z_final(5, npts)
    real(dp), intent(out) :: times_lost(npts)
    integer, intent(out) :: ierr(npts)

    integer :: i

    !$omp parallel do schedule(dynamic)
    do i = 1, npts
        call trace_orbit_soa_chunk(1, zstart(:, i:i), ntimstep, ntau, &
            dt_norm, ro0_in, atol, rtol, maxit, z_final(:, i:i), &
            times_lost(i:i), ierr(i:i))
    end do
    !$omp end parallel do
end subroutine trace_orbit_soa_omp1


subroutine trace_orbit_soa_chunk(npts, zstart, ntimstep, ntau, dt_norm, ro0_in, &
        atol, rtol, maxit, z_final, times_lost, ierr)
    use boozer_sub, only: vmec_to_boozer
    integer, intent(in) :: npts, ntimstep, ntau, maxit
    real(dp), intent(in) :: zstart(5, npts)
    real(dp), intent(in) :: dt_norm, ro0_in, atol, rtol
    real(dp), intent(out) :: z_final(5, npts)
    real(dp), intent(out) :: times_lost(npts)
    integer, intent(out) :: ierr(npts)

    real(dp) :: z_r(BATCH_PTS), z_th(BATCH_PTS), z_ph(BATCH_PTS), z_pphi(BATCH_PTS)
    real(dp) :: mu(BATCH_PTS), pabs(BATCH_PTS), vpar(BATCH_PTS)
    real(dp) :: Ath(BATCH_PTS), Aph(BATCH_PTS), dAth_dr(BATCH_PTS)
    real(dp) :: dAph_dr(BATCH_PTS), d2Aph_dr2(BATCH_PTS)
    real(dp) :: hth(BATCH_PTS), hph(BATCH_PTS)
    real(dp) :: dhth(3, BATCH_PTS), dhph(3, BATCH_PTS)
    real(dp) :: d2hth_dr2(BATCH_PTS), d2hph_dr2(BATCH_PTS)
    real(dp) :: Bmod(BATCH_PTS), dBmod(3, BATCH_PTS), d2Bmod(6, BATCH_PTS)
    logical :: escaped(BATCH_PTS)
    real(dp) :: ro0_norm, dt, theta_B, phi_B
    integer :: i, it

    if (npts <= 0) return

    ro0_norm = ro0_in / dsqrt(2.0d0)
    dt = dt_norm / dsqrt(2.0d0)

    do i = 1, npts
        z_r(i) = zstart(1, i)
        call vmec_to_boozer(zstart(1, i), mod(zstart(2, i), 6.283185307179586d0), &
            mod(zstart(3, i), 6.283185307179586d0), theta_B, phi_B)
        z_th(i) = mod(theta_B, 6.283185307179586d0)
        z_ph(i) = mod(phi_B, 6.283185307179586d0)
        pabs(i) = zstart(4, i)
    end do

    call eval_field_booz_many(npts, z_r, z_th, z_ph, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod)

    do i = 1, npts
        mu(i) = 0.5d0 * pabs(i)**2 * (1.0d0 - zstart(5, i)**2) / Bmod(i) * 2.0d0
        vpar(i) = pabs(i) * zstart(5, i) * dsqrt(2.0d0)
        z_pphi(i) = vpar(i) * hph(i) + Aph(i) / ro0_norm
    end do

    escaped(1:npts) = .false.
    ierr = 0
    times_lost = 0.0d0

    do it = 1, ntimstep
        call orbit_timestep_euler1_soa(npts, dt, ntau, ro0_norm, mu, &
            atol, rtol, maxit, z_r, z_th, z_ph, z_pphi, escaped, ierr)

        do i = 1, npts
            if (escaped(i) .and. times_lost(i) < 1.0d-30) then
                times_lost(i) = real(it * ntau, dp) * dt_norm
            end if
        end do

        if (all(escaped(1:npts))) exit
    end do

    call eval_field_booz_many(npts, z_r, z_th, z_ph, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod)

    do i = 1, npts
        z_final(1, i) = z_r(i)
        z_final(2, i) = z_th(i)
        z_final(3, i) = z_ph(i)
        vpar(i) = (z_pphi(i) - Aph(i) / ro0_norm) / hph(i)
        z_final(4, i) = dsqrt(mu(i) * Bmod(i) + 0.5d0 * vpar(i)**2)
        z_final(5, i) = vpar(i) / (z_final(4, i) * dsqrt(2.0d0))
    end do
end subroutine trace_orbit_soa_chunk

end module orbit_symplectic_soa
