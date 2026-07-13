module interface_crossing_ownership_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
contains
    subroutine stepped_field(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

        bmod = merge(2.0_dp, 1.0_dp, x(1) > 1.0_dp)
        sqrtg = 1.0_dp
        bder = 0.0_dp
        hcovar = [0.0_dp, 0.0_dp, 1.0_dp]
        hctrvr = hcovar
        hcurl = 0.0_dp
    end subroutine stepped_field
end module interface_crossing_ownership_fixture

program test_interface_crossing_ownership
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_sub, only: magfie
    use interface_crossing, only: apply_crossing, crossing_info_t, &
        CROSSING_LEVEL0, CROSS_CROSSING, CROSS_REFLECTION
    use interface_crossing_ownership_fixture, only: stepped_field

    implicit none

    real(dp) :: y_in(5), y_out(5)
    type(crossing_info_t) :: info

    magfie => stepped_field
    y_in = [1.0_dp, 0.2_dp, 0.3_dp, 1.0_dp, 1.0_dp]
    call apply_crossing(y_in, 1, 1, 3, CROSSING_LEVEL0, y_out, info)
    if (info%event_type /= CROSS_CROSSING) error stop 'outward crossing type'
    if (info%vol_to /= 2) error stop 'outward crossing ownership'
    if (y_out(1) /= 1.0_dp) error stop 'outward crossing moved off interface'

    call apply_crossing(y_in, 1, -1, 3, CROSSING_LEVEL0, y_out, info)
    if (info%event_type /= CROSS_CROSSING) error stop 'inward crossing type'
    if (info%vol_to /= 1) error stop 'inward crossing ownership'
    if (y_out(1) /= 1.0_dp) error stop 'inward crossing moved off interface'

    y_in(5) = 0.5_dp
    call apply_crossing(y_in, 1, 1, 3, CROSSING_LEVEL0, y_out, info)
    if (info%event_type /= CROSS_REFLECTION) error stop 'reflection type'
    if (info%vol_to /= 1) error stop 'reflection ownership'
    if (y_out(1) /= 1.0_dp) error stop 'reflection moved off interface'
end program test_interface_crossing_ownership
