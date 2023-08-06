!#include "assert.f90.inc"

program test_utility
use test_utils
use test_results
    double precision :: eps_flt = 0.1
    double precision :: eps_dbl = 1e-13
    
    logical :: test(6)
    
    test(1) = assert_int(1,1)                                                !Expect: pass
    test(2) = assert_int(1,2)                                                !Expect: fail
    test(3) = assert_flt(0.303, 0.202)                                 !Expect: fail
    test(4) = assert_flt_range(0.303, 0.300, eps_flt)    !Expect: pass
    test(5) = assert_flt(1e-12, 2e-12)                                 !Expect: pass (due to default eps = 1e-10)
    test(6) = assert_flt_range(1e-12, 2e-12, eps_dbl)   !Expect: fail
    
    call result_printer(test, 6, "Utility sanity test (expect o..oo.)")

end program test_utility
