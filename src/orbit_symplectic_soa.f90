module orbit_symplectic_soa
use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_boozer, only: eval_field_booz_many
use vector_potentail_mod, only: torflux

implicit none

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

    !$acc parallel loop
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
    !$acc end parallel loop
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

    !$acc parallel loop
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
    !$acc end parallel loop
end subroutine get_derivatives2_many


subroutine f_sympl_euler1_many(npts, dt, z_th, z_ph, pthold, ro0, mu, &
        x_r, x_pphi, fvec, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod, &
        pth, dpth, H, dH)
    integer, intent(in) :: npts
    real(dp), intent(in) :: dt, ro0
    real(dp), intent(in) :: mu(npts)
    real(dp), intent(in) :: z_th(npts), z_ph(npts), pthold(npts)
    real(dp), intent(in) :: x_r(npts), x_pphi(npts)
    real(dp), intent(out) :: fvec(2, npts)
    real(dp), intent(in) :: Ath(npts), Aph(npts)
    real(dp), intent(in) :: dAth_dr(npts), dAph_dr(npts), d2Aph_dr2(npts)
    real(dp), intent(in) :: hth(npts), hph(npts)
    real(dp), intent(in) :: dhth(3, npts), dhph(3, npts)
    real(dp), intent(in) :: d2hth_dr2(npts), d2hph_dr2(npts)
    real(dp), intent(in) :: Bmod(npts), dBmod(3, npts), d2Bmod(6, npts)
    real(dp), intent(out) :: pth(npts), dpth(4, npts)
    real(dp), intent(out) :: H(npts), dH(4, npts)

    real(dp), allocatable :: vpar(:), dvpar(:,:), d2vpar(:,:)
    real(dp), allocatable :: d2pth(:,:), d2H(:,:)
    integer :: i

    allocate(vpar(npts), dvpar(4, npts), d2vpar(10, npts))
    allocate(d2pth(10, npts), d2H(10, npts))

    call get_derivatives2_many(npts, ro0, mu, x_pphi, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod, &
        vpar, dvpar, d2vpar, pth, dpth, d2pth, H, dH, d2H)

    !$acc parallel loop
    do i = 1, npts
        fvec(1, i) = dpth(1, i) * (pth(i) - pthold(i)) &
            + dt * (dH(2, i) * dpth(1, i) - dH(1, i) * dpth(2, i))
        fvec(2, i) = dpth(1, i) * (x_pphi(i) - z_ph(i)) &
            + dt * (dH(3, i) * dpth(1, i) - dH(1, i) * dpth(3, i))
    end do
    !$acc end parallel loop
end subroutine f_sympl_euler1_many


subroutine jac_sympl_euler1_many(npts, dt, z_pphi, pthold, ro0, mu, &
        x_r, x_pphi, fjac, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod)
    integer, intent(in) :: npts
    real(dp), intent(in) :: dt, ro0
    real(dp), intent(in) :: mu(npts)
    real(dp), intent(in) :: z_pphi(npts), pthold(npts)
    real(dp), intent(in) :: x_r(npts), x_pphi(npts)
    real(dp), intent(out) :: fjac(2, 2, npts)
    real(dp), intent(in) :: Ath(npts), Aph(npts)
    real(dp), intent(in) :: dAth_dr(npts), dAph_dr(npts), d2Aph_dr2(npts)
    real(dp), intent(in) :: hth(npts), hph(npts)
    real(dp), intent(in) :: dhth(3, npts), dhph(3, npts)
    real(dp), intent(in) :: d2hth_dr2(npts), d2hph_dr2(npts)
    real(dp), intent(in) :: Bmod(npts), dBmod(3, npts), d2Bmod(6, npts)

    real(dp), allocatable :: vpar(:), dvpar(:,:), d2vpar(:,:)
    real(dp), allocatable :: pth(:), dpth(:,:), d2pth(:,:)
    real(dp), allocatable :: H(:), dH(:,:), d2H(:,:)
    integer :: i

    allocate(vpar(npts), dvpar(4, npts), d2vpar(10, npts))
    allocate(pth(npts), dpth(4, npts), d2pth(10, npts))
    allocate(H(npts), dH(4, npts), d2H(10, npts))

    call get_derivatives2_many(npts, ro0, mu, x_pphi, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod, &
        vpar, dvpar, d2vpar, pth, dpth, d2pth, H, dH, d2H)

    !$acc parallel loop
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
    !$acc end parallel loop
end subroutine jac_sympl_euler1_many


subroutine newton1_soa(npts, dt, ro0, mu, atol, rtol, maxit, &
        z_r, z_th, z_ph, z_pphi, pthold, x_r, x_pphi, converged)
    integer, intent(in) :: npts, maxit
    real(dp), intent(in) :: dt, ro0, atol, rtol
    real(dp), intent(in) :: mu(npts)
    real(dp), intent(in) :: z_r(npts), z_th(npts), z_ph(npts), z_pphi(npts)
    real(dp), intent(in) :: pthold(npts)
    real(dp), intent(inout) :: x_r(npts), x_pphi(npts)
    logical, intent(out) :: converged(npts)

    real(dp), allocatable :: Ath(:), Aph(:), dAth_dr(:), dAph_dr(:), d2Aph_dr2(:)
    real(dp), allocatable :: hth(:), hph(:), dhth(:,:), dhph(:,:)
    real(dp), allocatable :: d2hth_dr2(:), d2hph_dr2(:)
    real(dp), allocatable :: Bmod(:), dBmod(:,:), d2Bmod(:,:)
    real(dp), allocatable :: pth(:), dpth(:,:), H(:), dH(:,:)
    real(dp), allocatable :: fvec(:,:), fjac(:,:,:)
    real(dp), allocatable :: xlast_r(:), xlast_pphi(:)
    real(dp), allocatable :: tolref_pphi(:)
    real(dp) :: det, ijac11, ijac12, ijac21, ijac22
    real(dp) :: dx_r, dx_pphi
    integer :: kit, i

    allocate(Ath(npts), Aph(npts), dAth_dr(npts), dAph_dr(npts), d2Aph_dr2(npts))
    allocate(hth(npts), hph(npts), dhth(3, npts), dhph(3, npts))
    allocate(d2hth_dr2(npts), d2hph_dr2(npts))
    allocate(Bmod(npts), dBmod(3, npts), d2Bmod(6, npts))
    allocate(pth(npts), dpth(4, npts), H(npts), dH(4, npts))
    allocate(fvec(2, npts), fjac(2, 2, npts))
    allocate(xlast_r(npts), xlast_pphi(npts))
    allocate(tolref_pphi(npts))

    converged = .false.
    tolref_pphi = abs(10.0d0 * torflux / ro0)

    do kit = 1, maxit
        call eval_field_booz_many(npts, x_r, z_th, z_ph, &
            Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
            hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
            Bmod, dBmod, d2Bmod)

        call f_sympl_euler1_many(npts, dt, z_th, z_ph, pthold, ro0, mu, &
            x_r, x_pphi, fvec, &
            Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
            hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
            Bmod, dBmod, d2Bmod, &
            pth, dpth, H, dH)

        call jac_sympl_euler1_many(npts, dt, z_pphi, pthold, ro0, mu, &
            x_r, x_pphi, fjac, &
            Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
            hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
            Bmod, dBmod, d2Bmod)

        xlast_r = x_r
        xlast_pphi = x_pphi

        do i = 1, npts
            if (converged(i)) cycle

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

            if (abs(fvec(1, i)) < atol .and. abs(fvec(2, i)) < atol) then
                converged(i) = .true.
            else if (abs(x_r(i) - xlast_r(i)) < rtol .and. &
                     abs(x_pphi(i) - xlast_pphi(i)) < rtol * tolref_pphi(i)) then
                converged(i) = .true.
            end if
        end do

        if (all(converged)) exit
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

    real(dp), allocatable :: pthold(:), x_r(:), x_pphi(:)
    real(dp), allocatable :: Ath(:), Aph(:), dAth_dr(:), dAph_dr(:), d2Aph_dr2(:)
    real(dp), allocatable :: hth(:), hph(:), dhth(:,:), dhph(:,:)
    real(dp), allocatable :: d2hth_dr2(:), d2hph_dr2(:)
    real(dp), allocatable :: Bmod(:), dBmod(:,:), d2Bmod(:,:)
    real(dp), allocatable :: vpar(:), dvpar(:,:), pth(:), dpth(:,:), H(:), dH(:,:)
    logical, allocatable :: converged(:)
    integer :: ktau, i

    allocate(pthold(npts), x_r(npts), x_pphi(npts))
    allocate(Ath(npts), Aph(npts), dAth_dr(npts), dAph_dr(npts), d2Aph_dr2(npts))
    allocate(hth(npts), hph(npts), dhth(3, npts), dhph(3, npts))
    allocate(d2hth_dr2(npts), d2hph_dr2(npts))
    allocate(Bmod(npts), dBmod(3, npts), d2Bmod(6, npts))
    allocate(vpar(npts), dvpar(4, npts), pth(npts), dpth(4, npts), H(npts), dH(4, npts))
    allocate(converged(npts))

    escaped = .false.
    ierr = 0

    call eval_field_booz_many(npts, z_r, z_th, z_ph, &
        Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
        hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
        Bmod, dBmod, d2Bmod)

    call get_derivatives_many(npts, ro0, mu, z_pphi, &
        Ath, Aph, dAth_dr, dAph_dr, hth, hph, dhth, dhph, Bmod, dBmod, &
        vpar, dvpar, pth, dpth, H, dH)

    do ktau = 1, ntau
        pthold = pth

        x_r = z_r
        x_pphi = z_pphi

        call newton1_soa(npts, dt, ro0, mu, atol, rtol, maxit, &
            z_r, z_th, z_ph, z_pphi, pthold, x_r, x_pphi, converged)

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
        end do

        call eval_field_booz_many(npts, z_r, z_th, z_ph, &
            Ath, Aph, dAth_dr, dAph_dr, d2Aph_dr2, &
            hth, hph, dhth, dhph, d2hth_dr2, d2hph_dr2, &
            Bmod, dBmod, d2Bmod)

        call get_derivatives_many(npts, ro0, mu, z_pphi, &
            Ath, Aph, dAth_dr, dAph_dr, hth, hph, dhth, dhph, Bmod, dBmod, &
            vpar, dvpar, pth, dpth, H, dH)

        do i = 1, npts
            if (escaped(i)) cycle
            z_th(i) = z_th(i) + dt * dH(1, i) / dpth(1, i)
            z_ph(i) = z_ph(i) + dt * (vpar(i) - dH(1, i) / dpth(1, i) * hth(i)) / hph(i)
        end do
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

    real(dp), allocatable :: z_r(:), z_th(:), z_ph(:), z_pphi(:)
    real(dp), allocatable :: mu(:), pabs(:), vpar(:)
    real(dp), allocatable :: Ath(:), Aph(:), dAth_dr(:), dAph_dr(:), d2Aph_dr2(:)
    real(dp), allocatable :: hth(:), hph(:), dhth(:,:), dhph(:,:)
    real(dp), allocatable :: d2hth_dr2(:), d2hph_dr2(:)
    real(dp), allocatable :: Bmod(:), dBmod(:,:), d2Bmod(:,:)
    logical, allocatable :: escaped(:)
    real(dp) :: ro0_norm, dt, theta_B, phi_B
    integer :: i, it

    allocate(z_r(npts), z_th(npts), z_ph(npts), z_pphi(npts))
    allocate(mu(npts), pabs(npts), vpar(npts))
    allocate(Ath(npts), Aph(npts), dAth_dr(npts), dAph_dr(npts), d2Aph_dr2(npts))
    allocate(hth(npts), hph(npts), dhth(3, npts), dhph(3, npts))
    allocate(d2hth_dr2(npts), d2hph_dr2(npts))
    allocate(Bmod(npts), dBmod(3, npts), d2Bmod(6, npts))
    allocate(escaped(npts))

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

    escaped = .false.
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

        if (all(escaped)) exit
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
end subroutine trace_orbit_soa

end module orbit_symplectic_soa
