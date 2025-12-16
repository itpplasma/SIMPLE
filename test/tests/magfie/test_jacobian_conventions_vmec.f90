program test_jacobian_conventions_vmec
    !> Verify Jacobian conventions for VMEC reference coordinates.
    !>
    !> Checks:
    !> 1) signed Jacobian from covariant_basis triple product equals magfie_vmec sqrtg
    !> 2) abs(magfie_vmec sqrtg) equals sqrt(det g) from metric_tensor

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use magfie_sub, only: magfie_vmec
    use libneo_coordinates, only: coordinate_system_t, make_vmec_coordinate_system

    implicit none

    class(coordinate_system_t), allocatable :: vmec_cs
    real(dp) :: dummy

    integer, parameter :: n_r = 5
    integer, parameter :: n_th = 8
    integer, parameter :: n_ph = 6

    real(dp), parameter :: tol_signed = 5.0e-12_dp
    real(dp), parameter :: tol_abs = 5.0e-12_dp

    real(dp) :: r, theta, phi, phi_period
    real(dp) :: x_s(3)
    real(dp) :: bmod, sqrtg_vmec
    real(dp) :: bder(3), hcov(3), hctr(3), hcurl(3)
    real(dp) :: e_cov(3, 3)
    real(dp) :: g(3, 3), ginv(3, 3), sqrtg_metric
    real(dp) :: jac_signed
    real(dp) :: max_re_signed, max_re_abs
    integer :: i_r, i_th, i_ph
    integer :: n_failed

    n_failed = 0
    max_re_signed = 0.0_dp
    max_re_abs = 0.0_dp

    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call make_vmec_coordinate_system(vmec_cs)

    phi_period = twopi/real(nper, dp)

    do i_r = 1, n_r
        r = 0.25_dp + 0.6_dp*real(i_r - 1, dp)/real(n_r - 1, dp)
        do i_th = 1, n_th
            theta = 0.1_dp + (twopi - 0.2_dp)*real(i_th - 1, dp)/real(n_th - 1, dp)
            do i_ph = 1, n_ph
                phi = 0.1_dp + (phi_period - 0.2_dp)*real(i_ph - 1, dp)/ &
                    real(n_ph - 1, dp)

                x_s = [r**2, theta, phi]

                call magfie_vmec(x_s, bmod, sqrtg_vmec, bder, hcov, hctr, hcurl)

                call vmec_cs%covariant_basis(x_s, e_cov)
                jac_signed = signed_jacobian(e_cov)

                call vmec_cs%metric_tensor(x_s, g, ginv, sqrtg_metric)

                call update_max_scalar(jac_signed, sqrtg_vmec, max_re_signed)
                call update_max_scalar(abs(sqrtg_vmec), sqrtg_metric, max_re_abs)

                if (.not. approx_rel_or_abs(jac_signed, sqrtg_vmec, tol_signed, &
                                            1.0e-14_dp)) then
                    call report_fail('signed sqrtg', x_s, jac_signed, sqrtg_vmec, &
                                     tol_signed, n_failed)
                end if

                if (.not. approx_rel_or_abs(abs(sqrtg_vmec), sqrtg_metric, tol_abs, &
                                            1.0e-14_dp)) then
                    call report_fail('abs sqrtg', x_s, abs(sqrtg_vmec), sqrtg_metric, &
                                     tol_abs, n_failed)
                end if
            end do
        end do
    end do

    if (n_failed > 0) then
        print *, 'Max relative error signed sqrtg vs magfie_vmec:', max_re_signed
        print *, 'Max relative error abs(sqrtg) vs metric_tensor:', max_re_abs
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' Jacobian convention mismatch(es)'
        print *, '================================'
        error stop 1
    end if

    print *, 'Max relative error signed sqrtg vs magfie_vmec:', max_re_signed
    print *, 'Max relative error abs(sqrtg) vs metric_tensor:', max_re_abs
    print *, '================================'
    print *, 'PASSED: Jacobian conventions match for VMEC reference coordinates'
    print *, '================================'

contains

    pure logical function approx_rel_or_abs(a, b, tol_rel, tol_abs) result(ok)
        real(dp), intent(in) :: a, b, tol_rel, tol_abs

        if (abs(b) <= tol_abs) then
            ok = abs(a - b) <= tol_abs
        else
            ok = abs(a - b)/abs(b) <= tol_rel
        end if
    end function approx_rel_or_abs


    subroutine report_fail(name, x, a, b, tol, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3), a, b, tol
        integer, intent(inout) :: n_failed
        real(dp) :: rel

        rel = abs(a - b)/max(abs(b), 1.0e-30_dp)
        print *, 'FAIL: ', trim(name), ' at x=(s,th,ph)=', x
        print *, '  got=', a, ' expected=', b, ' rel=', rel, ' tol=', tol
        n_failed = n_failed + 1
    end subroutine report_fail


    subroutine update_max_scalar(a, b, max_rel)
        real(dp), intent(in) :: a, b
        real(dp), intent(inout) :: max_rel
        real(dp) :: denom, rel

        denom = max(abs(b), 1.0e-14_dp)
        rel = abs(a - b)/denom
        max_rel = max(max_rel, rel)
    end subroutine update_max_scalar


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

end program test_jacobian_conventions_vmec
