! Add anything wanted for test operations.
module test_utils
    implicit none
    save
    
    public assert_equal, check, check_close, check_close_vec, check_string
    interface assert_equal
        module procedure :: assert_int, assert_flt, assert_flt_range
        module procedure :: assert_dbl, assert_dbl_range
    end interface assert_equal

    contains
    
    
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                       Equality Templates
#define T int
#define TT integer
#include "templates/assert.f90.inc"

#define T flt
#define TT real
#include "templates/assert_flt.f90.inc"

#define T flt_range
#define TT real
#include "templates/assert_range.f90.inc"

#define T dbl
#define TT double precision
#include "templates/assert_flt.f90.inc"

#define T dbl_range
#define TT double precision
#include "templates/assert_range.f90.inc"


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                       Further Templates
subroutine check(condition, message, errors)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer, intent(inout) :: errors

    if (.not. condition) then
        print *, "ERROR:", trim(message)
        errors = errors + 1
    end if
end subroutine check

subroutine check_close(actual, expected, tolerance, message, errors)
    double precision, intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    integer, intent(inout) :: errors

    if (abs(actual - expected) > tolerance) then
        print *, "ERROR:", trim(message)
        print *, "  actual=", actual, " expected=", expected
        errors = errors + 1
    end if
end subroutine check_close

subroutine check_close_vec(actual, expected, tolerance, message, errors)
    double precision, intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: message
    integer, intent(inout) :: errors

    if (size(actual) /= size(expected)) then
        print *, "ERROR:", trim(message)
        print *, "  actual size=", size(actual), " expected size=", size(expected)
        errors = errors + 1
    elseif (maxval(abs(actual - expected)) > tolerance) then
        print *, "ERROR:", trim(message)
        print *, "  actual=", actual
        print *, "  expected=", expected
        errors = errors + 1
    end if
end subroutine check_close_vec

subroutine check_string(actual, expected, message, errors)
    character(len=*), intent(in) :: actual, expected, message
    integer, intent(inout) :: errors

    if (trim(actual) /= trim(expected)) then
        print *, "ERROR:", trim(message)
        print *, "  actual=", trim(actual), " expected=", trim(expected)
        errors = errors + 1
    end if
end subroutine check_string

end module test_utils
