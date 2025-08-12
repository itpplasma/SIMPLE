program test_orbit_symplectic_base_refactored
    use orbit_symplectic_base
    implicit none

    ! Test results tracking
    integer :: total_tests = 0, passed_tests = 0
    
    call run_all_tests()
    call print_test_summary()
    
contains

    subroutine run_all_tests()
        call test_init_rk_coefficients()
        call test_is_valid_method()
        call test_init_symplectic_integrator()
        call test_init_multistage_integrator()
    end subroutine run_all_tests

    ! Test: init_rk_coefficients for different methods
    subroutine test_init_rk_coefficients()
        real(dp), dimension(3,3) :: a, ahat
        real(dp), dimension(3) :: b, c
        real(dp), dimension(1,1) :: a1, ahat1
        real(dp), dimension(1) :: b1, c1
        logical :: test_passed
        
        total_tests = total_tests + 1
        test_passed = .true.
        
        ! Given: LOBATTO3 method with 3 stages
        ! When: initializing RK coefficients
        call init_rk_coefficients(LOBATTO3, 3, a, b, c, ahat)
        
        ! Then: coefficients should match expected Lobatto IIIA-IIIB values
        if (abs(a(1,1) - 0.0d0) > 1.0d-12) test_passed = .false.
        if (abs(b(1) - 0.16666666666666667d0) > 1.0d-12) test_passed = .false.
        if (abs(c(1) - 0.0d0) > 1.0d-12) test_passed = .false.
        if (abs(ahat(1,1) - 0.16666666666666667d0) > 1.0d-12) test_passed = .false.
        
        if (test_passed) then
            passed_tests = passed_tests + 1
            print *, "PASS: init_rk_coefficients - LOBATTO3 coefficients correct"
        else
            print *, "FAIL: init_rk_coefficients - LOBATTO3 coefficients incorrect"
        end if
        
        ! Test GAUSS1 method
        total_tests = total_tests + 1
        test_passed = .true.
        
        call init_rk_coefficients(GAUSS1, 1, a1, b1, c1)
        
        if (abs(a1(1,1) - 0.5d0) > 1.0d-12) test_passed = .false.
        if (abs(b1(1) - 1.0d0) > 1.0d-12) test_passed = .false.
        if (abs(c1(1) - 0.5d0) > 1.0d-12) test_passed = .false.
        
        if (test_passed) then
            passed_tests = passed_tests + 1
            print *, "PASS: init_rk_coefficients - GAUSS1 coefficients correct"
        else
            print *, "FAIL: init_rk_coefficients - GAUSS1 coefficients incorrect"
        end if
    end subroutine test_init_rk_coefficients

    ! Test: is_valid_method for different method-stage combinations
    subroutine test_is_valid_method()
        logical :: test_passed
        
        total_tests = total_tests + 1
        test_passed = .true.
        
        ! Given: different method-stage combinations
        ! When: checking validity
        ! Then: should return correct validation results
        
        if (.not. is_valid_method(GAUSS1, 1)) test_passed = .false.
        if (is_valid_method(GAUSS1, 2)) test_passed = .false.  ! Invalid combination
        if (.not. is_valid_method(LOBATTO3, 3)) test_passed = .false.
        if (is_valid_method(LOBATTO3, 2)) test_passed = .false.  ! Invalid combination
        if (.not. is_valid_method(RK45, 5)) test_passed = .false.  ! RK45 doesn't use stage count
        
        if (test_passed) then
            passed_tests = passed_tests + 1
            print *, "PASS: is_valid_method - validation logic correct"
        else
            print *, "FAIL: is_valid_method - validation logic incorrect"
        end if
    end subroutine test_is_valid_method

    ! Test: init_symplectic_integrator initialization
    subroutine test_init_symplectic_integrator()
        type(SymplecticIntegrator) :: si
        real(dp), dimension(4) :: z0
        real(dp) :: dt, atol, rtol
        logical :: test_passed
        
        total_tests = total_tests + 1
        test_passed = .true.
        
        ! Initialize test values
        z0 = [1.0d0, 0.5d0, 0.0d0, 2.0d0]
        dt = 0.01d0
        atol = 1.0d-12
        rtol = 1.0d-10
        
        ! Given: initial conditions and tolerances
        ! When: initializing symplectic integrator
        call init_symplectic_integrator(si, z0, dt, atol, rtol)
        
        ! Then: integrator should be properly initialized
        if (any(abs(si%z - z0) > 1.0d-15)) test_passed = .false.
        if (abs(si%dt - dt) > 1.0d-15) test_passed = .false.
        if (abs(si%atol - atol) > 1.0d-15) test_passed = .false.
        if (abs(si%rtol - rtol) > 1.0d-15) test_passed = .false.
        if (abs(si%pthold - 0.0d0) > 1.0d-15) test_passed = .false.
        if (si%ntau /= 1) test_passed = .false.
        if (abs(si%pabs - 1.0d0) > 1.0d-15) test_passed = .false.
        
        if (test_passed) then
            passed_tests = passed_tests + 1
            print *, "PASS: init_symplectic_integrator - initialization correct"
        else
            print *, "FAIL: init_symplectic_integrator - initialization incorrect"
        end if
    end subroutine test_init_symplectic_integrator

    ! Test: init_multistage_integrator initialization
    subroutine test_init_multistage_integrator()
        type(MultistageIntegrator) :: mi
        logical :: test_passed
        
        total_tests = total_tests + 1
        test_passed = .true.
        
        ! Given: valid method and stage count
        ! When: initializing multistage integrator
        call init_multistage_integrator(mi, 2, GAUSS2)
        
        ! Then: multistage integrator should be properly initialized
        if (mi%s /= 2) test_passed = .false.
        if (abs(mi%alpha(1) - 0.5d0) > 1.0d-15) test_passed = .false.
        if (abs(mi%beta(1) - 0.5d0) > 1.0d-15) test_passed = .false.
        
        if (test_passed) then
            passed_tests = passed_tests + 1
            print *, "PASS: init_multistage_integrator - initialization correct"
        else
            print *, "FAIL: init_multistage_integrator - initialization incorrect"
        end if
    end subroutine test_init_multistage_integrator

    subroutine print_test_summary()
        print *
        print *, "============================================"
        print *, "Test Summary for orbit_symplectic_base_refactored"
        print *, "============================================"
        print *, "Total tests: ", total_tests
        print *, "Passed: ", passed_tests
        print *, "Failed: ", total_tests - passed_tests
        if (passed_tests == total_tests) then
            print *, "Result: ALL TESTS PASSED"
        else
            print *, "Result: SOME TESTS FAILED"
        end if
        print *, "============================================"
    end subroutine print_test_summary

end program test_orbit_symplectic_base_refactored