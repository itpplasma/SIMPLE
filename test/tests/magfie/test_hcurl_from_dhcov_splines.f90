program test_hcurl_from_dhcov_splines
    !> Verify hcurl from derivatives of hcovar in reference coordinates.
    !>
    !> Uses analytic dhcov from splined_field_t%evaluate_with_der, and compares
    !> to a 5-point finite-difference oracle based on splined_field_t%evaluate.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use libneo_coordinates, only: coordinate_system_t

    implicit none

    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: f
    real(dp) :: dummy

    integer, parameter :: n_points = 240
    real(dp), parameter :: h_rho = 1.0e-4_dp
    real(dp), parameter :: h_ang = 1.0e-5_dp

    real(dp), parameter :: tol_rel = 1.0e-6_dp
    real(dp), parameter :: tol_abs = 1.0e-12_dp

    real(dp) :: phi_period
    real(dp) :: max_rel(3), max_abs(3)
    integer :: n_failed, i

    n_failed = 0
    max_rel = 0.0_dp
    max_abs = 0.0_dp

    call init_vmec('wout_ncsx.nc', 5, 5, 5, dummy)
    call init_reference_coordinates('wout_ncsx.nc')

    call create_coils_field('coils.simple', raw_coils)
    call create_splined_field(raw_coils, ref_coords, f)

    phi_period = twopi/real(nper, dp)

    do i = 1, n_points
        call check_point(i, phi_period, f, max_abs, max_rel, n_failed)
    end do

    call print_summary(max_abs, max_rel)

    if (n_failed > 0) then
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' hcurl mismatch(es)'
        print *, '================================'
        error stop 1
    end if

    print *, '================================'
    print *, 'PASSED: hcurl from dhcov matches FD oracle'
    print *, '================================'

contains

    subroutine check_point(i, phi_period, f, max_abs, max_rel, n_failed)
        integer, intent(in) :: i
        real(dp), intent(in) :: phi_period
        type(splined_field_t), intent(in) :: f
        real(dp), intent(inout) :: max_abs(3), max_rel(3)
        integer, intent(inout) :: n_failed

        real(dp) :: x(3)
        real(dp) :: Acov(3), hcov(3), Bmod
        real(dp) :: dAcov(3, 3), dhcov(3, 3), dBmod(3)
        real(dp) :: dhcov_fd(3, 3)
        real(dp) :: hcurl_new(3), hcurl_fd(3)
        real(dp) :: sqrtg

        x(1) = 0.25_dp + 0.55_dp*real(mod(13*i, 997), dp)/997.0_dp
        x(2) = 0.1_dp + (twopi - 0.2_dp)*real(mod(29*i, 991), dp)/991.0_dp
        x(3) = 0.1_dp + (phi_period - 0.2_dp)*real(mod(31*i, 983), dp)/983.0_dp

        call f%evaluate_with_der(x, Acov, hcov, Bmod, dAcov, dhcov, dBmod)

        call sqrtg_signed_scaled(ref_coords, x, sqrtg)

        call hcov_derivatives_fd5(f, x, phi_period, dhcov_fd)

        call compute_hcurl(sqrtg, dhcov, hcurl_new)
        call compute_hcurl(sqrtg, dhcov_fd, hcurl_fd)

        call update_and_check_vec('hcurl', x, hcurl_new, hcurl_fd, max_abs, &
                                  max_rel, n_failed)
    end subroutine check_point


    subroutine hcov_derivatives_fd5(f, x, phi_period, dh)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3), phi_period
        real(dp), intent(out) :: dh(3, 3)

        real(dp) :: hp1(3), hm1(3), hp2(3), hm2(3)
        integer :: k
        real(dp) :: step

        do k = 1, 3
            call eval_h_only(f, x, k, +1, phi_period, hp1)
            call eval_h_only(f, x, k, -1, phi_period, hm1)
            call eval_h_only(f, x, k, +2, phi_period, hp2)
            call eval_h_only(f, x, k, -2, phi_period, hm2)

            if (k == 1) then
                step = h_rho
            else
                step = h_ang
            end if

            dh(k, :) = (-hp2 + 8.0_dp*hp1 - 8.0_dp*hm1 + hm2)/(12.0_dp*step)
        end do
    end subroutine hcov_derivatives_fd5


    subroutine eval_h_only(f, x, k, m, phi_period, hcov)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: k, m
        real(dp), intent(in) :: phi_period
        real(dp), intent(out) :: hcov(3)

        real(dp) :: x_shift(3)
        real(dp) :: Acov(3), Bmod

        x_shift = x
        select case (k)
        case (1)
            x_shift(1) = x(1) + real(m, dp)*h_rho
        case (2)
            x_shift(2) = wrap_angle(x(2) + real(m, dp)*h_ang, twopi)
        case (3)
            x_shift(3) = wrap_angle(x(3) + real(m, dp)*h_ang, phi_period)
        case default
            error stop 'eval_h_only: invalid coordinate index'
        end select

        call f%evaluate(x_shift, Acov, hcov, Bmod)
    end subroutine eval_h_only


    subroutine compute_hcurl(sqrtg, dh, hcurl)
        real(dp), intent(in) :: sqrtg
        real(dp), intent(in) :: dh(3, 3)
        real(dp), intent(out) :: hcurl(3)

        if (abs(sqrtg) <= 0.0_dp) error stop 'compute_hcurl: sqrtg must be nonzero'

        hcurl(1) = (dh(2, 3) - dh(3, 2))/sqrtg
        hcurl(2) = (dh(3, 1) - dh(1, 3))/sqrtg
        hcurl(3) = (dh(1, 2) - dh(2, 1))/sqrtg
    end subroutine compute_hcurl


    subroutine sqrtg_signed_scaled(cs, x_scaled, sqrtg)
        class(coordinate_system_t), intent(in) :: cs
        real(dp), intent(in) :: x_scaled(3)
        real(dp), intent(out) :: sqrtg

        real(dp) :: x_ref(3)
        real(dp) :: e_cov(3, 3)
        real(dp) :: jac_signed, J

        ! Default scaling for VMEC reference coordinates in SIMPLE splines is
        ! x_scaled(1) = rho = sqrt(s), so s = rho^2 and ds/drho = 2*rho.
        x_ref = [x_scaled(1)**2, x_scaled(2), x_scaled(3)]
        J = 2.0_dp*x_scaled(1)

        call cs%covariant_basis(x_ref, e_cov)

        jac_signed = signed_jacobian(e_cov)
        sqrtg = jac_signed*J
    end subroutine sqrtg_signed_scaled


    pure function signed_jacobian(e_cov) result(jac)
        real(dp), intent(in) :: e_cov(3, 3)
        real(dp) :: jac
        real(dp) :: c(3)

        c = cross_product(e_cov(:, 2), e_cov(:, 3))
        jac = dot_product(e_cov(:, 1), c)
    end function signed_jacobian


    pure function cross_product(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product


    subroutine update_and_check_vec(name, x, a, b, max_abs, max_rel, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3)
        real(dp), intent(in) :: a(3), b(3)
        real(dp), intent(inout) :: max_abs(3), max_rel(3)
        integer, intent(inout) :: n_failed

        integer :: k
        real(dp) :: abs_err, denom, rel

        do k = 1, 3
            abs_err = abs(a(k) - b(k))
            denom = max(abs(b(k)), tol_abs)
            rel = abs_err/denom

            max_abs(k) = max(max_abs(k), abs_err)
            max_rel(k) = max(max_rel(k), rel)

            if (abs_err > tol_abs + tol_rel*abs(b(k))) then
                print *, 'FAIL: ', trim(name), '(', k, ') at x=(rho,th,ph)=', x
                print *, '  got=', a(k), ' expected=', b(k), ' abs=', abs_err, &
                         ' rel=', rel
                n_failed = n_failed + 1
            end if
        end do
    end subroutine update_and_check_vec


    subroutine print_summary(max_abs, max_rel)
        real(dp), intent(in) :: max_abs(3), max_rel(3)

        print *, 'Max abs errors in hcurl components:'
        print '(a,3(1x,es12.4))', '  hcurl_abs:', max_abs
        print *, 'Max rel errors in hcurl components:'
        print '(a,3(1x,es12.4))', '  hcurl_rel:', max_rel
    end subroutine print_summary


    pure real(dp) function wrap_angle(angle, period) result(wrapped)
        real(dp), intent(in) :: angle, period

        wrapped = modulo(angle, period)
        if (wrapped < 0.0_dp) wrapped = wrapped + period
    end function wrap_angle

end program test_hcurl_from_dhcov_splines
