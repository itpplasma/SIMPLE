program test_get_canonical_coordinates_refactored
! Unit tests for refactored get_canonical_coordinates module
! Tests follow behavior-driven design with Given-When-Then structure

use get_can_sub, only: dp, print_progress, interpolate_polynomial_backward, &
                       compute_derivative_factors, compute_derivative_factors_2nd
implicit none

integer :: n_tests_passed, n_tests_failed
real(dp), parameter :: tol = 1.0e-12_dp

n_tests_passed = 0
n_tests_failed = 0

call test_print_progress_utility()
call test_polynomial_interpolation()
call test_derivative_factors()
call test_second_derivative_factors()

print *, '============================================'
print *, 'Test Summary:'
print *, 'Passed:', n_tests_passed
print *, 'Failed:', n_tests_failed
print *, '============================================'

if (n_tests_failed > 0) then
  error stop 'Some tests failed!'
end if

contains

  !===========================================================================
  ! Test: Progress utility prints correctly
  !===========================================================================
  subroutine test_print_progress_utility()
    ! Given: Progress tracking parameters
    ! When: We call print_progress with different values
    ! Then: The function should execute without errors
    
    print *, 'Testing print_progress utility...'
    
    ! Test various progress values
    call print_progress('Test: ', 1, 10)
    call print_progress('Test: ', 5, 10)
    call print_progress('Test: ', 10, 10)
    
    print *, '  PASSED: Progress utility executes correctly'
    n_tests_passed = n_tests_passed + 1
    
  end subroutine test_print_progress_utility

  !===========================================================================
  ! Test: Polynomial interpolation backward method
  !===========================================================================
  subroutine test_polynomial_interpolation()
    ! Given: A simple polynomial data array
    ! When: We interpolate at different points
    ! Then: The results should match expected values
    
    real(dp), parameter :: data(4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    real(dp) :: result, expected
    real(dp) :: ds
    logical :: test_passed
    
    print *, 'Testing polynomial interpolation...'
    
    test_passed = .true.
    
    ! Test interpolation at ds = 0 (should give first coefficient)
    ds = 0.0_dp
    call interpolate_polynomial_backward(data, 3, ds, result)
    expected = 1.0_dp  ! data(1) + 0*(data(2) + 0*(data(3) + 0*data(4)))
    if (abs(result - expected) > tol) then
      test_passed = .false.
      print *, '  Failed at ds=0: expected', expected, 'got', result
    end if
    
    ! Test interpolation at ds = 1
    ds = 1.0_dp
    call interpolate_polynomial_backward(data, 3, ds, result)
    expected = 10.0_dp  ! 1 + 1*(2 + 1*(3 + 1*4)) = 1 + 1*(2 + 1*7) = 1 + 9 = 10
    if (abs(result - expected) > tol) then
      test_passed = .false.
      print *, '  Failed at ds=1: expected', expected, 'got', result
    end if
    
    ! Test edge case with n_order = 0
    call interpolate_polynomial_backward(data, 0, ds, result)
    if (abs(result) > tol) then
      test_passed = .false.
      print *, '  Failed n_order=0 case: expected 0, got', result
    end if
    
    if (test_passed) then
      print *, '  PASSED: Polynomial interpolation works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Polynomial interpolation validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_polynomial_interpolation

  !===========================================================================
  ! Test: Derivative factors computation
  !===========================================================================
  subroutine test_derivative_factors()
    ! Given: Data array and derivative factors
    ! When: We compute function value and derivative
    ! Then: Results should be consistent with manual calculation
    
    real(dp), parameter :: data(4) = [2.0_dp, 3.0_dp, 1.0_dp, 4.0_dp]
    real(dp), parameter :: derf_factors(4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    real(dp) :: result, dresult, ds
    logical :: test_passed
    
    print *, 'Testing derivative factors computation...'
    
    test_passed = .true.
    ds = 0.5_dp
    
    call compute_derivative_factors(data, derf_factors, 3, ds, result, dresult)
    
    ! Verify that function doesn't crash and returns reasonable values
    if (result < 0.0_dp .or. result > 1000.0_dp) then
      test_passed = .false.
      print *, '  Function value out of reasonable range:', result
    end if
    
    if (abs(dresult) > 1000.0_dp) then
      test_passed = .false.
      print *, '  Derivative value out of reasonable range:', dresult
    end if
    
    ! Test edge case n_order = 0
    call compute_derivative_factors(data, derf_factors, 0, ds, result, dresult)
    if (abs(result) > tol .or. abs(dresult) > tol) then
      test_passed = .false.
      print *, '  Failed n_order=0: result=', result, 'dresult=', dresult
    end if
    
    if (test_passed) then
      print *, '  PASSED: Derivative factors computation works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Derivative factors computation validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_derivative_factors

  !===========================================================================
  ! Test: Second derivative factors computation
  !===========================================================================
  subroutine test_second_derivative_factors()
    ! Given: Data and two derivative factor arrays
    ! When: We compute function, first, and second derivatives
    ! Then: The computation should complete without errors
    
    real(dp), parameter :: data(4) = [1.5_dp, 2.5_dp, 0.5_dp, 3.5_dp]
    real(dp), parameter :: derf1(4) = [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
    real(dp), parameter :: derf2(4) = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp]
    real(dp) :: result, dresult, d2result, ds
    logical :: test_passed
    
    print *, 'Testing second derivative factors computation...'
    
    test_passed = .true.
    ds = 0.3_dp
    
    call compute_derivative_factors_2nd(data, derf1, derf2, 3, ds, result, dresult, d2result)
    
    ! Check that computation completes and returns finite values
    if (.not. (result > -1000.0_dp .and. result < 1000.0_dp)) then
      test_passed = .false.
      print *, '  Function value not finite or out of range:', result
    end if
    
    if (.not. (abs(dresult) < 1000.0_dp)) then
      test_passed = .false.
      print *, '  First derivative out of range:', dresult
    end if
    
    if (.not. (abs(d2result) < 1000.0_dp)) then
      test_passed = .false.
      print *, '  Second derivative out of range:', d2result
    end if
    
    ! Test edge case n_order = 0
    call compute_derivative_factors_2nd(data, derf1, derf2, 0, ds, result, dresult, d2result)
    if (abs(result) > tol .or. abs(dresult) > tol .or. abs(d2result) > tol) then
      test_passed = .false.
      print *, '  Failed n_order=0 case'
    end if
    
    ! Test consistency: if we set derf2 = 0, d2result should be ~0
    call compute_derivative_factors_2nd(data, derf1, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 2, ds, result, dresult, d2result)
    if (abs(d2result) > tol * 10.0_dp) then
      test_passed = .false.
      print *, '  Consistency test failed: expected ~0 for d2result with zero derf2, got', d2result
    end if
    
    if (test_passed) then
      print *, '  PASSED: Second derivative factors computation works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Second derivative factors computation validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_second_derivative_factors

end program test_get_canonical_coordinates_refactored