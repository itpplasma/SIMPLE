program test_boozer_start_period
    !> A start point in any field period must keep its toroidal sector when
    !> mapped to Boozer integrator coordinates, and map back to itself.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use boozer_sub, only: get_boozer_coordinates
    use field_can_boozer, only: ref_to_integ_boozer, integ_to_ref_boozer
    use new_vmec_stuff_mod, only: nper

    implicit none

    real(dp), parameter :: twopi = 2.0_dp*3.1415926535897932384626433832795_dp
    real(dp), parameter :: tol = 1.0e-8_dp
    real(dp) :: fper, xref(3), xinteg(3), xback(3), shift, err
    integer :: k, j, n_failed

    n_failed = 0
    call init_vmec('wout.nc', 5, 5, 5, fper)
    call get_boozer_coordinates

    do k = 0, nper - 1
        do j = 1, 3
            xref = [0.15_dp + 0.3_dp*real(j - 1, dp), 0.4_dp + 1.7_dp*real(j, dp), &
                    0.2_dp*fper*real(j, dp) + real(k, dp)*fper]
            call ref_to_integ_boozer(xref, xinteg)
            call integ_to_ref_boozer(xinteg, xback)
            shift = angle_distance(xinteg(3), xref(3))
            err = max(abs(xback(1) - xref(1)), angle_distance(xback(2), xref(2)), &
                      angle_distance(xback(3), xref(3)))
            if (shift > 0.5_dp*fper .or. err > tol) then
                print *, 'FAILED: period', k, ' xref', xref, ' xinteg', xinteg, &
                    ' xback', xback
                n_failed = n_failed + 1
            end if
        end do
    end do

    if (n_failed > 0) then
        print *, 'FAILED: ', n_failed, ' start point(s) changed field period'
        error stop 1
    end if
    print *, 'PASSED: Boozer start mapping keeps the field period of', nper, &
        ' periods'

contains

    pure real(dp) function angle_distance(a, b)
        real(dp), intent(in) :: a, b
        angle_distance = abs(modulo(a - b + 0.5_dp*twopi, twopi) - 0.5_dp*twopi)
    end function angle_distance

end program test_boozer_start_period
