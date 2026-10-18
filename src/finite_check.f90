module finite_check
    !> Finiteness test and quiet NaN without the IEEE intrinsic modules.
    !>
    !> gfortran saves and restores the floating-point environment on entry to and
    !> exit from every procedure whose scope accesses an IEEE intrinsic module,
    !> including through a procedure-local USE of any module that imports one.
    !> In the orbit-tracing loop this cost up to 17% of the run time. The
    !> exponent-bit test below is exact for IEEE binary64 and is unaffected by
    !> -ffast-math, which may fold x /= x comparisons.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    implicit none
    private
    public :: is_finite, quiet_nan

    integer(int64), parameter :: EXPONENT_MASK = 2047_int64
    integer(int64), parameter :: QUIET_NAN_BITS = 9221120237041090560_int64

contains

    elemental logical function is_finite(x)
        real(dp), intent(in) :: x

        is_finite = iand(shiftr(transfer(x, 0_int64), 52), EXPONENT_MASK) &
            /= EXPONENT_MASK
    end function is_finite

    pure real(dp) function quiet_nan()
        quiet_nan = transfer(QUIET_NAN_BITS, 1.0_dp)
    end function quiet_nan

end module finite_check
