program test_hctr_from_metric_vmec
    !> Verify hctrvr = ginv @ hcovar in VMEC reference coordinates.

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

    real(dp), parameter :: tol_rel = 1.0e-12_dp
    real(dp), parameter :: tol_abs = 1.0e-13_dp

    real(dp) :: r, theta, phi, phi_period
    real(dp) :: x_s(3)
    real(dp) :: bmod, sqrtg
    real(dp) :: bder(3), hcov(3), hctr_ref(3), hcurl(3)
    real(dp) :: hctr_new(3)
    real(dp) :: g(3, 3), ginv(3, 3), sqrtg_metric
    real(dp) :: max_abs(3), max_rel(3)
    integer :: i_r, i_th, i_ph
    integer :: n_failed

    n_failed = 0
    max_abs = 0.0_dp
    max_rel = 0.0_dp

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

                call magfie_vmec(x_s, bmod, sqrtg, bder, hcov, hctr_ref, hcurl)
                call vmec_cs%metric_tensor(x_s, g, ginv, sqrtg_metric)

                hctr_new = matmul(ginv, hcov)

                call update_max_vec(hctr_new, hctr_ref, max_abs, max_rel)
                call check_vec('hctr', x_s, hctr_new, hctr_ref, tol_rel, tol_abs, &
                               n_failed)
            end do
        end do
    end do

    if (n_failed > 0) then
        print *, 'Max abs errors in hctr components:', max_abs
        print *, 'Max rel errors in hctr components:', max_rel
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' hctrvr mismatch(es)'
        print *, '================================'
        error stop 1
    end if

    print *, 'Max abs errors in hctr components:', max_abs
    print *, 'Max rel errors in hctr components:', max_rel
    print *, '================================'
    print *, 'PASSED: hctrvr matches ginv @ hcovar (VMEC ref coords)'
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


    subroutine report_fail(name, x, a, b, tol_rel, tol_abs, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3), a, b, tol_rel, tol_abs
        integer, intent(inout) :: n_failed
        real(dp) :: rel

        rel = abs(a - b)/max(abs(b), 1.0e-30_dp)
        print *, 'FAIL: ', trim(name), ' at x=(s,th,ph)=', x
        print *, '  got=', a, ' expected=', b, ' rel=', rel, ' tol_rel=', tol_rel, &
                 ' tol_abs=', tol_abs
        n_failed = n_failed + 1
    end subroutine report_fail


    subroutine check_vec(name, x, a, b, tol_rel, tol_abs, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3)
        real(dp), intent(in) :: a(3), b(3), tol_rel, tol_abs
        integer, intent(inout) :: n_failed
        integer :: i

        do i = 1, 3
            if (.not. approx_rel_or_abs(a(i), b(i), tol_rel, tol_abs)) then
                call report_fail(trim(name)//'('//trim(int_to_str(i))//')', x, &
                                 a(i), b(i), tol_rel, tol_abs, n_failed)
            end if
        end do
    end subroutine check_vec


    pure function int_to_str(i) result(s)
        integer, intent(in) :: i
        character(len=16) :: s

        write(s, '(i0)') i
    end function int_to_str


    subroutine update_max_scalar(a, b, max_abs, max_rel)
        real(dp), intent(in) :: a, b
        real(dp), intent(inout) :: max_abs, max_rel
        real(dp) :: denom, rel, abs_err

        abs_err = abs(a - b)
        denom = max(abs(b), 1.0e-14_dp)
        rel = abs_err/denom
        max_abs = max(max_abs, abs_err)
        max_rel = max(max_rel, rel)
    end subroutine update_max_scalar


    subroutine update_max_vec(a, b, max_abs, max_rel)
        real(dp), intent(in) :: a(3), b(3)
        real(dp), intent(inout) :: max_abs(3), max_rel(3)
        integer :: i

        do i = 1, 3
            call update_max_scalar(a(i), b(i), max_abs(i), max_rel(i))
        end do
    end subroutine update_max_vec

end program test_hctr_from_metric_vmec
