program test_finite_check
    !> is_finite and quiet_nan against the IEEE intrinsic module as oracle,
    !> over special values, the binary64 range edges and random bit patterns.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, &
        ieee_value, ieee_positive_inf, ieee_negative_inf, ieee_quiet_nan, &
        ieee_signaling_nan
    use finite_check, only: is_finite, quiet_nan
    implicit none

    real(dp) :: special(12), r
    integer(int64) :: bits
    integer :: i, nbad

    special = [0.0_dp, -0.0_dp, 1.0_dp, -1.0_dp, huge(1.0_dp), -huge(1.0_dp), &
               tiny(1.0_dp), tiny(1.0_dp)*epsilon(1.0_dp), &
               ieee_value(1.0_dp, ieee_positive_inf), &
               ieee_value(1.0_dp, ieee_negative_inf), &
               ieee_value(1.0_dp, ieee_quiet_nan), &
               ieee_value(1.0_dp, ieee_signaling_nan)]
    nbad = count(is_finite(special) .neqv. ieee_is_finite(special))

    call random_seed()
    do i = 1, 100000
        call random_number(r)
        bits = int(r*huge(1_int64), int64)
        if (mod(i, 2) == 0) bits = ior(bits, ishft(2047_int64, 52))
        r = transfer(bits, r)
        if (is_finite(r) .neqv. ieee_is_finite(r)) nbad = nbad + 1
    end do

    if (.not. ieee_is_nan(quiet_nan())) nbad = nbad + 1
    if (is_finite(quiet_nan())) nbad = nbad + 1

    if (nbad > 0) then
        write (error_unit, '(a,i0)') 'finite_check mismatches: ', nbad
        error stop 1
    end if
    print '(a)', 'test_finite_check: passed'
end program test_finite_check
