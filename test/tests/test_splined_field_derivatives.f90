program test_splined_field_derivatives
    !> Verify splined_field_t derivative API against finite differences of
    !> splined_field_t%evaluate in the same coordinates.

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
    integer :: n_failed

    integer, parameter :: n_points = 30
    real(dp), parameter :: tol_rel = 1.0e-7_dp
    real(dp), parameter :: tol_abs = 1.0e-9_dp
    real(dp), parameter :: h_rho = 1.0e-4_dp
    real(dp), parameter :: h_ang = 1.0e-5_dp

    real(dp) :: phi_period
    integer :: i

    n_failed = 0

    call init_vmec('wout_ncsx.nc', 5, 5, 5, dummy)
    call init_reference_coordinates('wout_ncsx.nc')

    call create_coils_field('coils.simple', raw_coils)
    call create_splined_field(raw_coils, ref_coords, splined)

    phi_period = twopi/real(nper, dp)

    do i = 1, n_points
        call check_point(i, phi_period, splined, n_failed)
    end do

    if (n_failed > 0) then
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' derivative mismatch(es) in splined_field_t'
        print *, '================================'
        error stop 1
    end if

    print *, '================================'
    print *, 'PASSED: splined_field_t derivatives match finite differences'
    print *, '================================'

contains

    subroutine check_point(i, phi_period, f, n_failed)
        integer, intent(in) :: i
        real(dp), intent(in) :: phi_period
        type(splined_field_t), intent(in) :: f
        integer, intent(inout) :: n_failed

        real(dp) :: x(3), x_plus(3), x_minus(3)
        real(dp) :: Acov(3), hcov(3), Bmod
        real(dp) :: dAcov(3, 3), dhcov(3, 3), dBmod(3)
        real(dp) :: Acov_p1(3), hcov_p1(3), Bmod_p1
        real(dp) :: Acov_m1(3), hcov_m1(3), Bmod_m1
        real(dp) :: Acov_p2(3), hcov_p2(3), Bmod_p2
        real(dp) :: Acov_m2(3), hcov_m2(3), Bmod_m2
        real(dp) :: dAcov_fd(3, 3), dhcov_fd(3, 3), dBmod_fd(3)

        integer :: k

        x(1) = 0.25_dp + 0.55_dp*real(mod(13*i, 97), dp)/97.0_dp
        x(2) = 0.1_dp + (twopi - 0.2_dp)*real(mod(29*i, 101), dp)/101.0_dp
        x(3) = 0.1_dp + (phi_period - 0.2_dp)*real(mod(31*i, 103), dp)/103.0_dp

        call f%evaluate_with_der(x, Acov, hcov, Bmod, dAcov, dhcov, dBmod)

        do k = 1, 3
            call eval_shifted(f, x, k, +1, phi_period, Acov_p1, hcov_p1, Bmod_p1)
            call eval_shifted(f, x, k, -1, phi_period, Acov_m1, hcov_m1, Bmod_m1)
            call eval_shifted(f, x, k, +2, phi_period, Acov_p2, hcov_p2, Bmod_p2)
            call eval_shifted(f, x, k, -2, phi_period, Acov_m2, hcov_m2, Bmod_m2)

            if (k == 1) then
                dAcov_fd(k, :) = (-Acov_p2 + 8.0_dp*Acov_p1 - 8.0_dp*Acov_m1 + &
                                  Acov_m2)/(12.0_dp*h_rho)
                dhcov_fd(k, :) = (-hcov_p2 + 8.0_dp*hcov_p1 - 8.0_dp*hcov_m1 + &
                                  hcov_m2)/(12.0_dp*h_rho)
                dBmod_fd(k) = (-Bmod_p2 + 8.0_dp*Bmod_p1 - 8.0_dp*Bmod_m1 + &
                               Bmod_m2)/(12.0_dp*h_rho)
            else
                dAcov_fd(k, :) = (-Acov_p2 + 8.0_dp*Acov_p1 - 8.0_dp*Acov_m1 + &
                                  Acov_m2)/(12.0_dp*h_ang)
                dhcov_fd(k, :) = (-hcov_p2 + 8.0_dp*hcov_p1 - 8.0_dp*hcov_m1 + &
                                  hcov_m2)/(12.0_dp*h_ang)
                dBmod_fd(k) = (-Bmod_p2 + 8.0_dp*Bmod_p1 - 8.0_dp*Bmod_m1 + &
                               Bmod_m2)/(12.0_dp*h_ang)
            end if
        end do

        call check_tensor('dAcov', x, dAcov, dAcov_fd, n_failed)
        call check_tensor('dhcov', x, dhcov, dhcov_fd, n_failed)
        call check_vec('dBmod', x, dBmod, dBmod_fd, n_failed)
    end subroutine check_point


    subroutine check_tensor(name, x, a, b, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3)
        real(dp), intent(in) :: a(3, 3), b(3, 3)
        integer, intent(inout) :: n_failed

        integer :: i, j

        do j = 1, 3
            do i = 1, 3
                if (.not. approx_rel_or_abs(a(i, j), b(i, j), tol_rel, tol_abs)) then
                    call report_fail(trim(name)//'('//trim(int_to_str(i))//','// &
                                     trim(int_to_str(j))//')', x, a(i, j), b(i, j), &
                                     tol_rel, n_failed)
                end if
            end do
        end do
    end subroutine check_tensor


    subroutine check_vec(name, x, a, b, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3)
        real(dp), intent(in) :: a(3), b(3)
        integer, intent(inout) :: n_failed

        integer :: i

        do i = 1, 3
            if (.not. approx_rel_or_abs(a(i), b(i), tol_rel, tol_abs)) then
                call report_fail(trim(name)//'('//trim(int_to_str(i))//')', x, a(i), &
                                 b(i), tol_rel, n_failed)
            end if
        end do
    end subroutine check_vec


    pure logical function approx_rel_or_abs(a, b, tol_rel, tol_abs) result(ok)
        real(dp), intent(in) :: a, b, tol_rel, tol_abs

        if (abs(b) <= tol_abs) then
            ok = abs(a - b) <= tol_abs
        else
            ok = abs(a - b)/abs(b) <= tol_rel
        end if
    end function approx_rel_or_abs


    subroutine report_fail(name, x, a, b, tol_rel, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3), a, b, tol_rel
        integer, intent(inout) :: n_failed
        real(dp) :: rel

        rel = abs(a - b)/max(abs(b), tol_abs)
        print *, 'FAIL: ', trim(name), ' at x=(rho,th,ph)=', x
        print *, '  got=', a, ' expected=', b, ' rel=', rel, ' tol=', tol_rel
        n_failed = n_failed + 1
    end subroutine report_fail


    pure real(dp) function wrap_angle(angle, period) result(wrapped)
        real(dp), intent(in) :: angle, period

        wrapped = modulo(angle, period)
        if (wrapped < 0.0_dp) wrapped = wrapped + period
    end function wrap_angle


    pure function int_to_str(i) result(s)
        integer, intent(in) :: i
        character(len=16) :: s

        write(s, '(i0)') i
    end function int_to_str


    subroutine eval_shifted(f, x, k, m, phi_period, Acov, hcov, Bmod)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: k, m
        real(dp), intent(in) :: phi_period
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod

        real(dp) :: x_shift(3)

        x_shift = x
        select case (k)
        case (1)
            x_shift(1) = x(1) + real(m, dp)*h_rho
        case (2)
            x_shift(2) = wrap_angle(x(2) + real(m, dp)*h_ang, twopi)
        case (3)
            x_shift(3) = wrap_angle(x(3) + real(m, dp)*h_ang, phi_period)
        case default
            error stop 'eval_shifted: invalid coordinate index'
        end select

        call f%evaluate(x_shift, Acov, hcov, Bmod)
    end subroutine eval_shifted

end program test_splined_field_derivatives
