program test_boozer_spline_derivatives
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord
    use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
    use new_vmec_stuff_mod, only: nper

    implicit none

    integer, parameter :: n_points = 24
    real(dp), parameter :: twopi = 2.0_dp*3.1415926535897932384626433832795_dp
    real(dp), parameter :: tol_rel = 2.0e-6_dp
    real(dp), parameter :: tol2_rel = 5.0e-3_dp
    real(dp), parameter :: tol_abs = 1.0e-10_dp
    real(dp), parameter :: h1_s = 2.0e-6_dp
    real(dp), parameter :: h1_th = 2.0e-6_dp
    real(dp), parameter :: h1_ph = 2.0e-6_dp
    real(dp), parameter :: h2_s = 5.0e-5_dp
    real(dp), parameter :: h2_th = 5.0e-5_dp
    real(dp), parameter :: h2_ph = 5.0e-5_dp

    real(dp) :: fper
    real(dp) :: phi_period
    integer :: i
    integer :: n_failed

    n_failed = 0

    use_B_r = .true.
    use_del_tp_B = .false.

    call init_vmec('wout.nc', 5, 5, 5, fper)
    call get_boozer_coordinates

    phi_period = twopi/real(nper, dp)

    do i = 1, n_points
        call check_point(i, phi_period, n_failed)
    end do

    if (n_failed > 0) then
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' derivative mismatch(es)'
        print *, '================================'
        error stop 1
    end if

    print *, '================================'
    print *, 'PASSED: Boozer Bmod/Br derivatives match FD oracle'
    print *, '================================'

contains

    subroutine check_point(i, phi_period, n_failed)
        integer, intent(in) :: i
        real(dp), intent(in) :: phi_period
        integer, intent(inout) :: n_failed

        real(dp) :: s, th, ph
        real(dp) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, d2A_phi_ds2, d3A_phi_ds3
        real(dp) :: Bth, dBth, d2Bth
        real(dp) :: Bph, dBph, d2Bph
        real(dp) :: Bmod, Br
        real(dp) :: dBmod(3), dBr(3)
        real(dp) :: d2Bmod(6), d2Br(6)

        real(dp) :: dBmod_fd(3), dBr_fd(3)
        real(dp) :: d2Bmod_fd(3), d2Br_fd(3)

        s = 0.15_dp + 0.7_dp*real(mod(31*i, 997), dp)/997.0_dp
        th = 0.1_dp + (twopi - 0.2_dp)*real(mod(37*i, 991), dp)/991.0_dp
        ph = 0.1_dp + (phi_period - 0.2_dp)*real(mod(41*i, 983), dp)/983.0_dp

        call splint_boozer_coord(s, th, ph, A_theta, A_phi, dA_theta_ds, dA_phi_ds, &
                                 d2A_phi_ds2, d3A_phi_ds3, Bth, dBth, d2Bth, Bph, &
                                 dBph, d2Bph, Bmod, dBmod, d2Bmod, Br, dBr, d2Br)

        call fd_first_derivatives(s, th, ph, phi_period, Bmod, Br, dBmod_fd, dBr_fd)
        call fd_second_diag(s, th, ph, phi_period, Bmod, Br, d2Bmod_fd, d2Br_fd)

        call check_vec('dBmod', [dBmod(1), dBmod(2), dBmod(3)], dBmod_fd, &
                       tol_rel, s, th, ph, n_failed)
        call check_vec('dBr  ', [dBr(1), dBr(2), dBr(3)], dBr_fd, &
                       tol_rel, s, th, ph, n_failed)

        call check_vec('d2Bmod(diag)', [d2Bmod(1), d2Bmod(4), d2Bmod(6)], d2Bmod_fd, &
                       tol2_rel, s, th, ph, n_failed)
        call check_vec('d2Br  (diag)', [d2Br(1), d2Br(4), d2Br(6)], d2Br_fd, &
                       tol2_rel, s, th, ph, n_failed)
    end subroutine check_point


    subroutine fd_first_derivatives(s, th, ph, phi_period, Bmod0, Br0, dBmod, dBr)
        real(dp), intent(in) :: s, th, ph, phi_period, Bmod0, Br0
        real(dp), intent(out) :: dBmod(3), dBr(3)

        real(dp) :: Bp, Bm
        real(dp) :: Brp, Brm

        call eval_B_only(s + h1_s, th, ph, Bp, Brp)
        call eval_B_only(s - h1_s, th, ph, Bm, Brm)
        dBmod(1) = (Bp - Bm)/(2.0_dp*h1_s)
        dBr(1) = (Brp - Brm)/(2.0_dp*h1_s)

        call eval_B_only(s, wrap_angle(th + h1_th, twopi), ph, Bp, Brp)
        call eval_B_only(s, wrap_angle(th - h1_th, twopi), ph, Bm, Brm)
        dBmod(2) = (Bp - Bm)/(2.0_dp*h1_th)
        dBr(2) = (Brp - Brm)/(2.0_dp*h1_th)

        call eval_B_only(s, th, wrap_angle(ph + h1_ph, phi_period), Bp, Brp)
        call eval_B_only(s, th, wrap_angle(ph - h1_ph, phi_period), Bm, Brm)
        dBmod(3) = (Bp - Bm)/(2.0_dp*h1_ph)
        dBr(3) = (Brp - Brm)/(2.0_dp*h1_ph)

        if (.not. is_finite(Bmod0) .or. .not. is_finite(Br0)) then
            error stop 'fd_first_derivatives: non-finite baseline values'
        end if
    end subroutine fd_first_derivatives


    subroutine fd_second_diag(s, th, ph, phi_period, Bmod0, Br0, d2Bmod, d2Br)
        real(dp), intent(in) :: s, th, ph, phi_period, Bmod0, Br0
        real(dp), intent(out) :: d2Bmod(3), d2Br(3)

        real(dp) :: Bp1, Bm1, Bp2, Bm2
        real(dp) :: Brp1, Brm1, Brp2, Brm2

        call eval_B_only(s + h2_s, th, ph, Bp1, Brp1)
        call eval_B_only(s - h2_s, th, ph, Bm1, Brm1)
        call eval_B_only(s + 2.0_dp*h2_s, th, ph, Bp2, Brp2)
        call eval_B_only(s - 2.0_dp*h2_s, th, ph, Bm2, Brm2)
        d2Bmod(1) = (-Bp2 + 16.0_dp*Bp1 - 30.0_dp*Bmod0 + 16.0_dp*Bm1 - Bm2) / &
                    (12.0_dp*h2_s*h2_s)
        d2Br(1) = (-Brp2 + 16.0_dp*Brp1 - 30.0_dp*Br0 + 16.0_dp*Brm1 - Brm2) / &
                  (12.0_dp*h2_s*h2_s)

        call eval_B_only(s, wrap_angle(th + h2_th, twopi), ph, Bp1, Brp1)
        call eval_B_only(s, wrap_angle(th - h2_th, twopi), ph, Bm1, Brm1)
        call eval_B_only(s, wrap_angle(th + 2.0_dp*h2_th, twopi), ph, Bp2, Brp2)
        call eval_B_only(s, wrap_angle(th - 2.0_dp*h2_th, twopi), ph, Bm2, Brm2)
        d2Bmod(2) = (-Bp2 + 16.0_dp*Bp1 - 30.0_dp*Bmod0 + 16.0_dp*Bm1 - Bm2) / &
                    (12.0_dp*h2_th*h2_th)
        d2Br(2) = (-Brp2 + 16.0_dp*Brp1 - 30.0_dp*Br0 + 16.0_dp*Brm1 - Brm2) / &
                  (12.0_dp*h2_th*h2_th)

        call eval_B_only(s, th, wrap_angle(ph + h2_ph, phi_period), Bp1, Brp1)
        call eval_B_only(s, th, wrap_angle(ph - h2_ph, phi_period), Bm1, Brm1)
        call eval_B_only(s, th, wrap_angle(ph + 2.0_dp*h2_ph, phi_period), Bp2, Brp2)
        call eval_B_only(s, th, wrap_angle(ph - 2.0_dp*h2_ph, phi_period), Bm2, Brm2)
        d2Bmod(3) = (-Bp2 + 16.0_dp*Bp1 - 30.0_dp*Bmod0 + 16.0_dp*Bm1 - Bm2) / &
                    (12.0_dp*h2_ph*h2_ph)
        d2Br(3) = (-Brp2 + 16.0_dp*Brp1 - 30.0_dp*Br0 + 16.0_dp*Brm1 - Brm2) / &
                  (12.0_dp*h2_ph*h2_ph)
    end subroutine fd_second_diag


    subroutine eval_B_only(s, th, ph, Bmod, Br)
        real(dp), intent(in) :: s, th, ph
        real(dp), intent(out) :: Bmod, Br

        real(dp) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, d2A_phi_ds2, d3A_phi_ds3
        real(dp) :: Bth, dBth, d2Bth
        real(dp) :: Bph, dBph, d2Bph
        real(dp) :: dBmod(3), dBr(3)
        real(dp) :: d2Bmod(6), d2Br(6)

        call splint_boozer_coord(s, th, ph, A_theta, A_phi, dA_theta_ds, dA_phi_ds, &
                                 d2A_phi_ds2, d3A_phi_ds3, Bth, dBth, d2Bth, Bph, &
                                 dBph, d2Bph, Bmod, dBmod, d2Bmod, Br, dBr, d2Br)
    end subroutine eval_B_only


    subroutine check_vec(name, a, b, tol_rel, s, th, ph, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: a(3), b(3)
        real(dp), intent(in) :: tol_rel
        real(dp), intent(in) :: s, th, ph
        integer, intent(inout) :: n_failed

        integer :: k
        real(dp) :: rel, denom

        do k = 1, 3
            denom = max(abs(b(k)), tol_abs)
            rel = abs(a(k) - b(k))/denom
            if (rel > tol_rel) then
                print *, 'FAIL: ', trim(name), ' k=', k, ' at (s,th,ph)=', s, th, ph
                print *, '  got=', a(k), ' expected=', b(k), ' rel=', rel
                n_failed = n_failed + 1
            end if
        end do
    end subroutine check_vec


    pure real(dp) function wrap_angle(angle, period) result(wrapped)
        real(dp), intent(in) :: angle, period

        wrapped = modulo(angle, period)
        if (wrapped < 0.0_dp) wrapped = wrapped + period
    end function wrap_angle


    pure logical function is_finite(x) result(ok)
        real(dp), intent(in) :: x
        ok = (x == x) .and. (abs(x) < huge(x))
    end function is_finite

end program test_boozer_spline_derivatives
