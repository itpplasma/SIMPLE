! Add anything wanted for test operations.
module test_utils
    implicit none
    save
    
    public assert_equal
    interface assert_equal
        module procedure :: assert_int, assert_flt, assert_flt_range, assert_dbl, assert_dbl_range
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
end module test_utils
