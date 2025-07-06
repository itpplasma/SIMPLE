program test_array_utils
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use array_utils, only: init_derivative_factors
  implicit none
  
  integer :: i, errors
  
  errors = 0
  
  ! Test basic functionality
  call test_basic_values(errors)
  
  ! Test edge cases
  call test_edge_cases(errors)
  
  ! Test large arrays
  call test_large_arrays(errors)
  
  ! Test numerical accuracy
  call test_numerical_accuracy(errors)
  
  if (errors == 0) then
    print *, "All array_utils tests passed!"
  else
    print *, "ERROR: ", errors, " test(s) failed!"
    stop 1
  end if
  
contains

  subroutine test_basic_values(errors)
    integer, intent(inout) :: errors
    double precision :: derf1(10), derf2(10), derf3(10)
    double precision :: expected1(10), expected2(10), expected3(10)
    integer :: k
    
    print *, "Testing basic derivative factor values..."
    
    ! Initialize expected values
    do k = 1, 10
      expected1(k) = dble(k-1)
      expected2(k) = dble((k-1)*(k-2))
      expected3(k) = dble((k-1)*(k-2)*(k-3))
    end do
    
    ! Call the function
    call init_derivative_factors(10, derf1, derf2, derf3)
    
    ! Check results
    do k = 1, 10
      if (abs(derf1(k) - expected1(k)) > 1.0d-15) then
        print *, "ERROR: derf1(", k, ") = ", derf1(k), " expected ", expected1(k)
        errors = errors + 1
      end if
      if (abs(derf2(k) - expected2(k)) > 1.0d-15) then
        print *, "ERROR: derf2(", k, ") = ", derf2(k), " expected ", expected2(k)
        errors = errors + 1
      end if
      if (abs(derf3(k) - expected3(k)) > 1.0d-15) then
        print *, "ERROR: derf3(", k, ") = ", derf3(k), " expected ", expected3(k)
        errors = errors + 1
      end if
    end do
    
    if (errors == 0) then
      print *, "  Basic values test PASSED"
    end if
    
  end subroutine test_basic_values
  
  subroutine test_edge_cases(errors)
    integer, intent(inout) :: errors
    double precision :: derf1(5), derf2(5), derf3(5)
    
    print *, "Testing edge cases..."
    
    ! Test with small array
    call init_derivative_factors(5, derf1, derf2, derf3)
    
    ! Check specific edge values
    ! For k=1: derf1 = 0, derf2 = 0, derf3 = 0
    if (abs(derf1(1)) > 1.0d-15 .or. abs(derf2(1)) > 1.0d-15 .or. abs(derf3(1)) > 1.0d-15) then
      print *, "ERROR: k=1 should give all zeros"
      errors = errors + 1
    end if
    
    ! For k=2: derf1 = 1, derf2 = 0, derf3 = 0
    if (abs(derf1(2) - 1.0d0) > 1.0d-15 .or. abs(derf2(2)) > 1.0d-15 .or. abs(derf3(2)) > 1.0d-15) then
      print *, "ERROR: k=2 values incorrect"
      errors = errors + 1
    end if
    
    ! For k=3: derf1 = 2, derf2 = 2, derf3 = 0
    if (abs(derf1(3) - 2.0d0) > 1.0d-15 .or. abs(derf2(3) - 2.0d0) > 1.0d-15 .or. abs(derf3(3)) > 1.0d-15) then
      print *, "ERROR: k=3 values incorrect"
      errors = errors + 1
    end if
    
    ! For k=4: derf1 = 3, derf2 = 6, derf3 = 6
    if (abs(derf1(4) - 3.0d0) > 1.0d-15 .or. abs(derf2(4) - 6.0d0) > 1.0d-15 .or. abs(derf3(4) - 6.0d0) > 1.0d-15) then
      print *, "ERROR: k=4 values incorrect"
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Edge cases test PASSED"
    end if
    
  end subroutine test_edge_cases
  
  subroutine test_large_arrays(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n_large = 1000
    double precision :: derf1(n_large), derf2(n_large), derf3(n_large)
    integer :: k
    
    print *, "Testing large arrays..."
    
    ! Initialize large arrays
    call init_derivative_factors(n_large, derf1, derf2, derf3)
    
    ! Check some specific values
    ! For k=100: derf1 = 99, derf2 = 99*98 = 9702, derf3 = 99*98*97 = 941094
    if (abs(derf1(100) - 99.0d0) > 1.0d-15) then
      print *, "ERROR: derf1(100) incorrect"
      errors = errors + 1
    end if
    if (abs(derf2(100) - 9702.0d0) > 1.0d-15) then
      print *, "ERROR: derf2(100) incorrect"
      errors = errors + 1
    end if
    if (abs(derf3(100) - 941094.0d0) > 1.0d-15) then
      print *, "ERROR: derf3(100) incorrect"
      errors = errors + 1
    end if
    
    ! Check last value
    if (abs(derf1(n_large) - dble(n_large-1)) > 1.0d-15) then
      print *, "ERROR: derf1(", n_large, ") incorrect"
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Large arrays test PASSED"
    end if
    
  end subroutine test_large_arrays
  
  subroutine test_numerical_accuracy(errors)
    integer, intent(inout) :: errors
    double precision :: derf1(50), derf2(50), derf3(50)
    double precision :: factorial_ratio
    integer :: k
    
    print *, "Testing numerical accuracy..."
    
    call init_derivative_factors(50, derf1, derf2, derf3)
    
    ! Test relationships between arrays
    ! derf2(k) should equal derf1(k) * (k-2) for k >= 3
    do k = 3, 50
      factorial_ratio = derf2(k) / derf1(k)
      if (abs(factorial_ratio - dble(k-2)) > 1.0d-15) then
        print *, "ERROR: derf2/derf1 ratio incorrect at k=", k
        errors = errors + 1
      end if
    end do
    
    ! derf3(k) should equal derf2(k) * (k-3) for k >= 4
    do k = 4, 50
      factorial_ratio = derf3(k) / derf2(k)
      if (abs(factorial_ratio - dble(k-3)) > 1.0d-15) then
        print *, "ERROR: derf3/derf2 ratio incorrect at k=", k
        errors = errors + 1
      end if
    end do
    
    if (errors == 0) then
      print *, "  Numerical accuracy test PASSED"
    end if
    
  end subroutine test_numerical_accuracy

end program test_array_utils