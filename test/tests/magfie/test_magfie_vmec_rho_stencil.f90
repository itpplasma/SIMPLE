program test_magfie_vmec_rho_stencil
    !> Near the axis |B| is smooth in rho = sqrt(s), so 2 rho d(log|B|)/ds =
    !> d(log|B|)/drho tends to a finite, angle-dependent limit. A radial stencil
    !> in s with a fixed step much larger than s misses that limit (#754).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use magfie_sub, only: magfie_vmec

    implicit none

    real(dp), parameter :: s_pairs(2, 3) = reshape([1.0e-6_dp, 4.0e-6_dp, &
                                                    1.0e-5_dp, 4.0e-5_dp, &
                                                    1.0e-4_dp, 4.0e-4_dp], [2, 3])
    real(dp), parameter :: thetas(3) = [0.3_dp, 2.1_dp, 4.4_dp]
    real(dp) :: g(2), dummy
    integer :: i, j, k, n_failed

    n_failed = 0
    call init_vmec('wout.nc', 5, 5, 5, dummy)

    do i = 1, size(s_pairs, 2)
        do j = 1, size(thetas)
            do k = 1, 2
                g(k) = log_b_rho_derivative(s_pairs(k, i), thetas(j))
            end do
            ! Between rho and 2 rho the limit changes by O(rho) only.
            if (abs(g(1) - g(2)) > 0.05_dp*max(abs(g(1)), abs(g(2)))) then
                print *, 'd(log|B|)/drho not smooth at s =', s_pairs(:, i), &
                    ' theta =', thetas(j), ':', g
                n_failed = n_failed + 1
            end if
        end do
    end do

    if (n_failed > 0) error stop 1
    print *, 'test_magfie_vmec_rho_stencil PASSED'

contains

    real(dp) function log_b_rho_derivative(s, theta) result(g)
        real(dp), intent(in) :: s, theta
        real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call magfie_vmec([s, theta, 0.2_dp], bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        g = 2.0_dp*sqrt(s)*bder(1)
    end function log_b_rho_derivative

end program test_magfie_vmec_rho_stencil
