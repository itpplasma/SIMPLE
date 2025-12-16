program test_bder_from_splines
    !> Verify bder = d(log Bmod)/dx computed from spline derivatives.
    !> Reference is a high-order finite-difference oracle of Bmod from evaluate().

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use reference_coordinates, only: init_reference_coordinates, ref_coords

    implicit none

    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: splined

    real(dp) :: dummy
    real(dp) :: phi_period
    integer :: n_failed

    integer, parameter :: n_points = 240
    real(dp), parameter :: tol_rel = 5.0e-8_dp
    real(dp), parameter :: tol_abs = 1.0e-10_dp
    real(dp), parameter :: h_rho = 1.0e-4_dp
    real(dp), parameter :: h_ang = 1.0e-5_dp

    real(dp) :: max_rel(3)
    integer :: i

    n_failed = 0
    max_rel = 0.0_dp

    call init_vmec('wout_ncsx.nc', 5, 5, 5, dummy)
    call init_reference_coordinates('wout_ncsx.nc')

    call create_coils_field('coils.simple', raw_coils)
    call create_splined_field(raw_coils, ref_coords, splined)

    phi_period = twopi/real(nper, dp)

    do i = 1, n_points
        call check_point(i, phi_period, splined, max_rel, n_failed)
    end do

    if (n_failed > 0) then
        call print_summary(max_rel)
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' bder mismatch(es)'
        print *, '================================'
        error stop 1
    end if

    call print_summary(max_rel)
    print *, '================================'
    print *, 'PASSED: bder from splines matches FD oracle'
    print *, '================================'

contains

    subroutine check_point(i, phi_period, f, max_rel, n_failed)
        integer, intent(in) :: i
        real(dp), intent(in) :: phi_period
        type(splined_field_t), intent(in) :: f
        real(dp), intent(inout) :: max_rel(3)
        integer, intent(inout) :: n_failed

        real(dp) :: x(3)
        real(dp) :: Acov(3), hcov(3), Bmod
        real(dp) :: dAcov(3, 3), dhcov(3, 3), dBmod(3)
        real(dp) :: bder_spline(3), bder_fd(3)
        integer :: k

        x(1) = 0.25_dp + 0.55_dp*real(mod(13*i, 997), dp)/997.0_dp
        x(2) = 0.1_dp + (twopi - 0.2_dp)*real(mod(29*i, 991), dp)/991.0_dp
        x(3) = 0.1_dp + (phi_period - 0.2_dp)*real(mod(31*i, 983), dp)/983.0_dp

        call f%evaluate_with_der(x, Acov, hcov, Bmod, dAcov, dhcov, dBmod)
        bder_spline = dBmod/max(Bmod, 1.0e-30_dp)

        do k = 1, 3
            bder_fd(k) = dlogB_dx_fd5(f, x, k, phi_period)/max(Bmod, 1.0e-30_dp)
            call update_max_and_check(k, x, bder_spline(k), bder_fd(k), max_rel, &
                                      n_failed)
        end do
    end subroutine check_point


    real(dp) function dlogB_dx_fd5(f, x, k, phi_period) result(dBdx)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: k
        real(dp), intent(in) :: phi_period

        real(dp) :: Bp1, Bm1, Bp2, Bm2
        real(dp) :: step

        call eval_B_only(f, x, k, +1, phi_period, Bp1)
        call eval_B_only(f, x, k, -1, phi_period, Bm1)
        call eval_B_only(f, x, k, +2, phi_period, Bp2)
        call eval_B_only(f, x, k, -2, phi_period, Bm2)

        if (k == 1) then
            step = h_rho
        else
            step = h_ang
        end if

        dBdx = (-Bp2 + 8.0_dp*Bp1 - 8.0_dp*Bm1 + Bm2)/(12.0_dp*step)
    end function dlogB_dx_fd5


    subroutine eval_B_only(f, x, k, m, phi_period, Bmod)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: k, m
        real(dp), intent(in) :: phi_period
        real(dp), intent(out) :: Bmod

        real(dp) :: x_shift(3)
        real(dp) :: Acov(3), hcov(3)

        x_shift = x
        select case (k)
        case (1)
            x_shift(1) = x(1) + real(m, dp)*h_rho
        case (2)
            x_shift(2) = wrap_angle(x(2) + real(m, dp)*h_ang, twopi)
        case (3)
            x_shift(3) = wrap_angle(x(3) + real(m, dp)*h_ang, phi_period)
        case default
            error stop 'eval_B_only: invalid coordinate index'
        end select

        call f%evaluate(x_shift, Acov, hcov, Bmod)
    end subroutine eval_B_only


    subroutine update_max_and_check(k, x, a, b, max_rel, n_failed)
        integer, intent(in) :: k
        real(dp), intent(in) :: x(3), a, b
        real(dp), intent(inout) :: max_rel(3)
        integer, intent(inout) :: n_failed

        real(dp) :: rel, denom

        denom = max(abs(b), tol_abs)
        rel = abs(a - b)/denom
        max_rel(k) = max(max_rel(k), rel)

        if (rel > tol_rel) then
            print *, 'FAIL: bder(', k, ') at x=(rho,th,ph)=', x
            print *, '  got=', a, ' expected=', b, ' rel=', rel, ' tol=', tol_rel
            n_failed = n_failed + 1
        end if
    end subroutine update_max_and_check


    subroutine print_summary(max_rel)
        real(dp), intent(in) :: max_rel(3)

        print *, 'Max relative errors in bder components:'
        print '(a,3(1x,es12.4))', '  bder:', max_rel
    end subroutine print_summary


    pure real(dp) function wrap_angle(angle, period) result(wrapped)
        real(dp), intent(in) :: angle, period

        wrapped = modulo(angle, period)
        if (wrapped < 0.0_dp) wrapped = wrapped + period
    end function wrap_angle

end program test_bder_from_splines
