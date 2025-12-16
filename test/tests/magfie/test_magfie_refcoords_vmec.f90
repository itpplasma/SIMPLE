program test_magfie_refcoords_vmec
    !> Decoupling test: reproduce magfie_vmec outputs using only
    !> (1) magnetic_field_t%evaluate in reference coordinates and
    !> (2) coordinate_system_t%metric_tensor (plus a known scaling).
    !>
    !> This keeps the VMEC equilibrium and current magfie_vmec as oracle, but
    !> exercises the refcoords-style interface needed to generalize RK45.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use magfie_sub, only: magfie_vmec
    use field_vmec, only: vmec_field_t, create_vmec_field
    use libneo_coordinates, only: coordinate_system_t, make_vmec_coordinate_system

    implicit none

    type(vmec_field_t) :: field
    class(coordinate_system_t), allocatable :: vmec_cs
    real(dp) :: dummy

    integer, parameter :: n_r = 5
    integer, parameter :: n_th = 6
    integer, parameter :: n_ph = 4

    real(dp), parameter :: hs = 1.0e-3_dp
    real(dp), parameter :: tol_bmod = 2.0e-3_dp
    real(dp), parameter :: tol_sqrtg = 2.0e-3_dp
    real(dp), parameter :: tol_bder = 1.0e-2_dp
    real(dp), parameter :: tol_hcov = 1.0e-2_dp
    real(dp), parameter :: tol_hctr = 1.0e-2_dp
    real(dp), parameter :: tol_hcurl = 5.0e-2_dp

    real(dp) :: r, theta, phi, phi_period
    real(dp) :: x_r(3), x_s(3)
    real(dp) :: bmod_or, sqrtg_or
    real(dp) :: bder_or(3), hcov_or(3), hctr_or(3), hcurl_or(3)
    real(dp) :: bmod_fd, sqrtg_fd
    real(dp) :: bder_fd(3), hcov_fd(3), hctr_fd(3), hcurl_fd(3)
    real(dp) :: sqrtg_or_abs, sgn_sqrtg
    real(dp) :: g_s(3, 3), ginv_s(3, 3), ginv_r(3, 3), sqrtg_s
    real(dp) :: dh(3, 3), dBdx(3)
    real(dp) :: J, dr, ht, hp
    integer :: i_r, i_th, i_ph
    integer :: n_failed

    n_failed = 0
    phi_period = twopi/real(nper, dp)
    ht = hs*twopi
    hp = ht/5.0_dp

    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call create_vmec_field(field)
    call make_vmec_coordinate_system(vmec_cs)

    do i_r = 1, n_r
        r = 0.25_dp + 0.6_dp*real(i_r - 1, dp)/real(n_r - 1, dp)
        do i_th = 1, n_th
            theta = 0.1_dp + (twopi - 0.2_dp)*real(i_th - 1, dp)/real(n_th - 1, dp)
            do i_ph = 1, n_ph
                phi = 0.1_dp + (phi_period - 0.2_dp)*real(i_ph - 1, dp)/ &
                    real(n_ph - 1, dp)

                x_r = [r, theta, phi]
                x_s = [r**2, theta, phi]

                J = 2.0_dp*r
                dr = hs/J

                call vmec_cs%metric_tensor(x_s, g_s, ginv_s, sqrtg_s)
                call metric_inverse_scaled(J, ginv_s, ginv_r)
                sqrtg_fd = sqrtg_s*J

                call magfie_refcoords_fd(field, x_r, dr, ht, hp, phi_period, &
                                         bmod_fd, bder_fd, hcov_fd, hctr_fd, &
                                         hcurl_fd, dBdx, dh)

                call magfie_vmec(x_s, bmod_or, sqrtg_or, bder_or, hcov_or, hctr_or, &
                                 hcurl_or)
                call transform_magfie_s_to_r(J, bmod_or, sqrtg_or, bder_or, hcov_or, &
                                             hctr_or, hcurl_or)

                sgn_sqrtg = sign(1.0_dp, sqrtg_or)
                sqrtg_or_abs = abs(sqrtg_or)
                hcurl_or = hcurl_or*sgn_sqrtg

                if (.not. approx_rel_or_abs(bmod_fd, bmod_or, tol_bmod, &
                                            1.0e-12_dp)) then
                    call report_fail('bmod', x_r, bmod_fd, bmod_or, tol_bmod, n_failed)
                end if
                if (.not. approx_rel_or_abs(sqrtg_fd, sqrtg_or_abs, tol_sqrtg, &
                                            1.0e-10_dp)) then
                    call report_fail('sqrtg', x_r, sqrtg_fd, sqrtg_or_abs, &
                                     tol_sqrtg, n_failed)
                end if
                call check_vec('bder', x_r, bder_fd, bder_or, tol_bder, 1.0e-10_dp, &
                               n_failed)
                call check_vec('hcov', x_r, hcov_fd, hcov_or, tol_hcov, 1.0e-10_dp, &
                               n_failed)
                call check_vec('hctr', x_r, hctr_fd, hctr_or, tol_hctr, 1.0e-10_dp, &
                               n_failed)
                call check_vec('hcurl', x_r, hcurl_fd, hcurl_or, tol_hcurl, &
                               1.0e-10_dp, n_failed)
            end do
        end do
    end do

    if (n_failed > 0) then
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' mismatch(es) in refcoords magfie test'
        print *, '================================'
        error stop 1
    end if

    print *, '================================'
    print *, 'PASSED: magfie via refcoords interface matches magfie_vmec (VMEC ref)'
    print *, '================================'

contains

    subroutine magfie_refcoords_fd(field, x, dr, dth, dph, phi_period, bmod, bder, &
                                   hcov, hctr, hcurl, dBdx, dh)
        type(vmec_field_t), intent(in) :: field
        real(dp), intent(in) :: x(3), dr, dth, dph, phi_period
        real(dp), intent(out) :: bmod, bder(3), hcov(3), hctr(3), hcurl(3)
        real(dp), intent(out) :: dBdx(3), dh(3, 3)

        real(dp) :: B0, Bp, Bm
        real(dp) :: h0(3), hpv(3), hmv(3)
        real(dp) :: x_plus(3), x_minus(3)
        real(dp) :: Acov_tmp(3)
        integer :: i

        call field%evaluate(x, Acov_tmp, h0, B0)
        bmod = B0
        hcov = h0

        do i = 1, 3
            x_plus = x
            x_minus = x
            select case (i)
            case (1)
                x_plus(1) = x(1) + dr
                x_minus(1) = x(1) - dr
            case (2)
                x_plus(2) = wrap_angle(x(2) + dth, twopi)
                x_minus(2) = wrap_angle(x(2) - dth, twopi)
            case (3)
                x_plus(3) = wrap_angle(x(3) + dph, phi_period)
                x_minus(3) = wrap_angle(x(3) - dph, phi_period)
            case default
                error stop 'magfie_refcoords_fd: invalid coordinate index'
            end select

            call field%evaluate(x_plus, Acov_tmp, hpv, Bp)
            call field%evaluate(x_minus, Acov_tmp, hmv, Bm)

            select case (i)
            case (1)
                dBdx(i) = (Bp - Bm)/(2.0_dp*dr)
                dh(i, :) = (hpv - hmv)/(2.0_dp*dr)
            case (2)
                dBdx(i) = (Bp - Bm)/(2.0_dp*dth)
                dh(i, :) = (hpv - hmv)/(2.0_dp*dth)
            case (3)
                dBdx(i) = (Bp - Bm)/(2.0_dp*dph)
                dh(i, :) = (hpv - hmv)/(2.0_dp*dph)
            end select
        end do

        bder = dBdx/max(bmod, 1.0e-30_dp)

        call compute_hcurl(sqrtg_fd, dh, hcurl)

        hctr = matmul(ginv_r, hcov)
    end subroutine magfie_refcoords_fd


    subroutine metric_inverse_scaled(J, ginv_s, ginv_r)
        real(dp), intent(in) :: J
        real(dp), intent(in) :: ginv_s(3, 3)
        real(dp), intent(out) :: ginv_r(3, 3)

        ginv_r = ginv_s
        ginv_r(1, 1) = ginv_s(1, 1)/(J*J)
        ginv_r(1, 2) = ginv_s(1, 2)/J
        ginv_r(1, 3) = ginv_s(1, 3)/J
        ginv_r(2, 1) = ginv_s(2, 1)/J
        ginv_r(3, 1) = ginv_s(3, 1)/J
    end subroutine metric_inverse_scaled


    subroutine transform_magfie_s_to_r(J, bmod, sqrtg, bder, hcov, hctr, hcurl)
        real(dp), intent(in) :: J
        real(dp), intent(inout) :: bmod, sqrtg
        real(dp), intent(inout) :: bder(3), hcov(3), hctr(3), hcurl(3)

        sqrtg = sqrtg*J
        bder(1) = bder(1)*J
        hcov(1) = hcov(1)*J
        hctr(1) = hctr(1)/J
        hcurl(1) = hcurl(1)/J
    end subroutine transform_magfie_s_to_r


    pure real(dp) function wrap_angle(angle, period) result(wrapped)
        real(dp), intent(in) :: angle, period

        wrapped = modulo(angle, period)
        if (wrapped < 0.0_dp) wrapped = wrapped + period
    end function wrap_angle


    subroutine compute_hcurl(sqrtg, dh, hcurl)
        real(dp), intent(in) :: sqrtg
        real(dp), intent(in) :: dh(3, 3)
        real(dp), intent(out) :: hcurl(3)

        if (.not. (sqrtg > 0.0_dp)) error stop 'compute_hcurl: sqrtg must be positive'

        hcurl(1) = (dh(2, 3) - dh(3, 2))/sqrtg
        hcurl(2) = (dh(3, 1) - dh(1, 3))/sqrtg
        hcurl(3) = (dh(1, 2) - dh(2, 1))/sqrtg
    end subroutine compute_hcurl


    pure logical function approx_rel_or_abs(a, b, tol_rel, tol_abs) result(ok)
        real(dp), intent(in) :: a, b, tol_rel, tol_abs

        if (abs(b) <= tol_abs) then
            ok = abs(a - b) <= tol_abs
        else
            ok = abs(a - b)/abs(b) <= tol_rel
        end if
    end function approx_rel_or_abs


    subroutine check_vec(name, x, a, b, tol_rel, tol_abs, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3)
        real(dp), intent(in) :: a(3), b(3), tol_rel, tol_abs
        integer, intent(inout) :: n_failed
        integer :: i

        do i = 1, 3
            if (.not. approx_rel_or_abs(a(i), b(i), tol_rel, tol_abs)) then
                call report_fail(trim(name)//'('//trim(int_to_str(i))//')', x, &
                                 a(i), b(i), tol_rel, n_failed)
            end if
        end do
    end subroutine check_vec


    subroutine report_fail(name, x, a, b, tol, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3), a, b, tol
        integer, intent(inout) :: n_failed
        real(dp) :: rel

        rel = abs(a - b)/max(abs(b), 1.0e-30_dp)
        print *, 'FAIL: ', trim(name), ' at x=(r,th,ph)=', x
        print *, '  got=', a, ' expected=', b, ' rel=', rel, ' tol=', tol
        n_failed = n_failed + 1
    end subroutine report_fail


    pure function int_to_str(i) result(s)
        integer, intent(in) :: i
        character(len=16) :: s

        write(s, '(i0)') i
    end function int_to_str

end program test_magfie_refcoords_vmec
