program test_magfie_geoflux

    use, intrinsic :: iso_fortran_env, only : dp => real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use simple, only : init_vmec
    use magfie_sub, only : init_magfie, VMEC, magfie
    use new_vmec_stuff_mod, only : nper
    use util, only : twopi

    implicit none

    real(dp) :: fper
    real(dp) :: x(3)
    real(dp) :: bmod, sqrtg
    real(dp) :: bder(3), hcov(3), hctrvr(3), hcurl(3)
    integer :: i

    call init_vmec('EQDSK_I.geqdsk', 5, 5, 5, fper)
    if (nper /= 1) then
        error stop 'test_magfie_geoflux: nper should be 1 for GEQDSK'
    end if

    call init_magfie(VMEC)

    call random_seed()

    do i = 1, 128
        call random_number(x)
        x(1) = min(0.999_dp, x(1))
        x(2) = (x(2) - 0.5_dp) * twopi
        x(3) = x(3) * twopi
        call magfie(x, bmod, sqrtg, bder, hcov, hctrvr, hcurl)

        if (bmod <= 0.0_dp) then
            error stop 'test_magfie_geoflux: Bmod not positive'
        end if
        if (.not. ieee_is_finite(sqrtg)) then
            error stop 'test_magfie_geoflux: sqrtg not finite'
        end if
        if (.not. all(ieee_is_finite(hcov))) then
            error stop 'test_magfie_geoflux: hcov not finite'
        end if
        if (.not. all(ieee_is_finite(hctrvr))) then
            error stop 'test_magfie_geoflux: hctrvr not finite'
        end if
        if (.not. all(ieee_is_finite(hcurl))) then
            error stop 'test_magfie_geoflux: hcurl not finite'
        end if
    end do

end program test_magfie_geoflux
