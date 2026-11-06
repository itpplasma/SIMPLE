program test_magfie_vmec_axis_stencil
    !> Within hs = 1e-3 of the axis the radial finite differences of magfie_vmec
    !> must not sample s < 0. All outputs stay finite there, and the radial
    !> derivative of |B| is continuous at s = hs.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use simple, only: init_vmec
    use magfie_sub, only: magfie_vmec

    implicit none

    real(dp), parameter :: hs = 1.0e-3_dp
    real(dp), parameter :: s_near(3) = [9.99e-4_dp, 3.0e-4_dp, 2.0e-5_dp]
    real(dp), parameter :: thetas(3) = [0.3_dp, 2.1_dp, 4.4_dp]
    real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: bder_below, bder_above, dummy
    integer :: i, j, n_failed

    n_failed = 0
    call init_vmec('wout.nc', 5, 5, 5, dummy)

    do i = 1, size(s_near)
        do j = 1, size(thetas)
            call magfie_vmec([s_near(i), thetas(j), 0.2_dp], bmod, sqrtg, bder, &
                             hcovar, hctrvr, hcurl)
            if (.not. (ieee_is_finite(bmod) .and. ieee_is_finite(sqrtg) .and. &
                     all(ieee_is_finite(bder)) .and. all(ieee_is_finite(hcovar)) .and. &
                     all(ieee_is_finite(hctrvr)) .and. all(ieee_is_finite(hcurl)))) then
                print *, 'non-finite magfie_vmec output at s =', s_near(i), &
                    ' theta =', thetas(j)
                n_failed = n_failed + 1
            end if
        end do
    end do

    do j = 1, size(thetas)
        call magfie_vmec([hs*(1.0_dp - 1.0e-9_dp), thetas(j), 0.2_dp], bmod, sqrtg, &
                         bder, hcovar, hctrvr, hcurl)
        bder_below = bder(1)
        call magfie_vmec([hs, thetas(j), 0.2_dp], bmod, sqrtg, bder, hcovar, &
                         hctrvr, hcurl)
        bder_above = bder(1)
        if (abs(bder_below - bder_above) > 1.0e-2_dp*abs(bder_above)) then
            print *, 'd|B|/ds jumps at s = hs: one-sided', bder_below, &
                ' central', bder_above
            n_failed = n_failed + 1
        end if
    end do

    if (n_failed > 0) error stop 1
    print *, 'test_magfie_vmec_axis_stencil PASSED'
end program test_magfie_vmec_axis_stencil
