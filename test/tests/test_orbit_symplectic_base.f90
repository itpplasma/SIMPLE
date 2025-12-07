program test_orbit_symplectic_base
  use orbit_symplectic_base
  implicit none
  
  integer :: errors
  
  errors = 0
  
  ! Test integration method constants
  call test_integration_constants(errors)
  
  ! Test Runge-Kutta Gauss coefficients
  call test_rk_gauss_coefficients(errors)
  
  ! Test Runge-Kutta Lobatto coefficients  
  call test_rk_lobatto_coefficients(errors)
  
  ! Test symplectic_integrator_t type initialization
  call test_symplectic_integrator_type(errors)
  
  if (errors == 0) then
    print *, "All orbit_symplectic_base module tests passed!"
  else
    print *, "ERROR: ", errors, " test(s) failed!"
    stop 1
  end if
  
contains

  subroutine test_integration_constants(errors)
    integer, intent(inout) :: errors
    integer :: methods(7)
    integer :: i, j
    logical :: all_distinct
    
    print *, "Testing integration method constants..."
    
    ! Given: The module defines constants for different integration methods
    ! When: We use these constants
    ! Then: They should be distinct and within expected range
    
    ! Test that constants are distinct (no two methods have same value)
    
    methods = [RK45, EXPL_IMPL_EULER, IMPL_EXPL_EULER, MIDPOINT, GAUSS1, GAUSS2, LOBATTO3]
    
    all_distinct = .true.
    do i = 1, size(methods)
      do j = i+1, size(methods)
        if (methods(i) == methods(j)) then
          print *, "ERROR: Integration method constants not distinct"
          print *, "Method", i, "and", j, "have same value:", methods(i)
          errors = errors + 1
          all_distinct = .false.
          exit
        end if
      end do
      if (.not. all_distinct) exit
    end do
    
    ! Test that S_MAX is large enough for all methods
    if (S_MAX < maxval(methods)) then
      print *, "ERROR: S_MAX too small for defined methods"
      print *, "S_MAX:", S_MAX, "Max method value:", maxval(methods)
      errors = errors + 1
    end if
    
    ! Test that all constants are non-negative (sensible range)
    do i = 1, size(methods)
      if (methods(i) < 0) then
        print *, "ERROR: Integration method constant is negative:", methods(i)
        errors = errors + 1
      end if
    end do
    
    if (errors == 0) then
      print *, "  Integration constants test PASSED"
    end if
    
  end subroutine test_integration_constants
  
  subroutine test_rk_gauss_coefficients(errors)
    integer, intent(inout) :: errors
    
    print *, "Testing Runge-Kutta Gauss coefficients..."
    
    ! Test 1-stage Gauss method (order 2)
    call test_gauss_n1(errors)
    
    ! Test 2-stage Gauss method (order 4)
    call test_gauss_n2(errors)
    
    ! Test 3-stage Gauss method (order 6)
    call test_gauss_n3(errors)
    
    ! Test 4-stage Gauss method (order 8)
    call test_gauss_n4(errors)
    
    ! Test unsupported stage count
    call test_gauss_unsupported(errors)
    
    if (errors == 0) then
      print *, "  RK Gauss coefficients test PASSED"
    end if
    
  end subroutine test_rk_gauss_coefficients
  
  subroutine test_gauss_n1(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 1
    real(dp) :: a(n,n), b(n), c(n)
    real(dp), parameter :: tol = 1.0d-14
    
    ! Given: A 1-stage Gauss method
    ! When: We compute the coefficients
    ! Then: They should match the known 1-stage Gauss values
    
    call coeff_rk_gauss(n, a, b, c)
    
    ! Check a coefficients
    if (abs(a(1,1) - 0.5d0) > tol) then
      print *, "ERROR: 1-stage Gauss a(1,1) incorrect"
      print *, "Expected: 0.5, Got:", a(1,1)
      errors = errors + 1
    end if
    
    ! Check b coefficients
    if (abs(b(1) - 1.0d0) > tol) then
      print *, "ERROR: 1-stage Gauss b(1) incorrect"
      print *, "Expected: 1.0, Got:", b(1)
      errors = errors + 1
    end if
    
    ! Check c coefficients
    if (abs(c(1) - 0.5d0) > tol) then
      print *, "ERROR: 1-stage Gauss c(1) incorrect"
      print *, "Expected: 0.5, Got:", c(1)
      errors = errors + 1
    end if
    
  end subroutine test_gauss_n1
  
  subroutine test_gauss_n2(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 2
    real(dp) :: a(n,n), b(n), c(n)
    real(dp), parameter :: tol = 1.0d-14
    
    ! Given: A 2-stage Gauss method
    ! When: We compute the coefficients
    ! Then: They should satisfy Gauss method properties
    
    call coeff_rk_gauss(n, a, b, c)
    
    ! Check symmetry properties
    if (abs(a(1,1) - a(2,2)) > tol) then
      print *, "ERROR: 2-stage Gauss should have a(1,1) = a(2,2)"
      errors = errors + 1
    end if
    
    ! Check that b coefficients sum to 1
    if (abs(sum(b) - 1.0d0) > tol) then
      print *, "ERROR: 2-stage Gauss b coefficients should sum to 1"
      print *, "Sum:", sum(b)
      errors = errors + 1
    end if
    
    ! Check that b coefficients are symmetric
    if (abs(b(1) - b(2)) > tol) then
      print *, "ERROR: 2-stage Gauss b coefficients should be symmetric"
      errors = errors + 1
    end if
    
    ! Check c coefficient symmetry
    if (abs(c(1) + c(2) - 1.0d0) > tol) then
      print *, "ERROR: 2-stage Gauss c coefficients should sum to 1"
      errors = errors + 1
    end if
    
  end subroutine test_gauss_n2
  
  subroutine test_gauss_n3(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 3
    real(dp) :: a(n,n), b(n), c(n)
    real(dp), parameter :: tol = 1.0d-12  ! Slightly looser tolerance for 3-stage
    
    call coeff_rk_gauss(n, a, b, c)
    
    ! Check symmetry of diagonal elements
    if (abs(a(1,1) - a(3,3)) > tol) then
      print *, "ERROR: 3-stage Gauss should have a(1,1) = a(3,3)"
      errors = errors + 1
    end if
    
    ! Check that b coefficients sum to 1
    if (abs(sum(b) - 1.0d0) > tol) then
      print *, "ERROR: 3-stage Gauss b coefficients should sum to 1"
      print *, "Sum:", sum(b)
      errors = errors + 1
    end if
    
    ! Check symmetry of b coefficients
    if (abs(b(1) - b(3)) > tol) then
      print *, "ERROR: 3-stage Gauss b(1) should equal b(3)"
      errors = errors + 1
    end if
    
    ! Check that c(2) = 0.5 for symmetric methods
    if (abs(c(2) - 0.5d0) > tol) then
      print *, "ERROR: 3-stage Gauss c(2) should be 0.5"
      print *, "Got:", c(2)
      errors = errors + 1
    end if
    
  end subroutine test_gauss_n3
  
  subroutine test_gauss_n4(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 4
    real(dp) :: a(n,n), b(n), c(n)
    real(dp), parameter :: tol = 1.0d-12
    
    call coeff_rk_gauss(n, a, b, c)
    
    ! Check that b coefficients sum to 1
    if (abs(sum(b) - 1.0d0) > tol) then
      print *, "ERROR: 4-stage Gauss b coefficients should sum to 1"
      print *, "Sum:", sum(b)
      errors = errors + 1
    end if
    
    ! Check symmetry of b coefficients
    if (abs(b(1) - b(4)) > tol .or. abs(b(2) - b(3)) > tol) then
      print *, "ERROR: 4-stage Gauss b coefficients should be symmetric"
      errors = errors + 1
    end if
    
    ! Check that c coefficients are in [0,1]
    if (any(c < 0.0d0) .or. any(c > 1.0d0)) then
      print *, "ERROR: 4-stage Gauss c coefficients should be in [0,1]"
      errors = errors + 1
    end if
    
  end subroutine test_gauss_n4
  
  subroutine test_gauss_unsupported(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 5  ! Unsupported stage count
    real(dp) :: a(n,n), b(n), c(n)
    real(dp), parameter :: tol = 1.0d-14
    
    ! Given: An unsupported stage count
    ! When: We call coeff_rk_gauss
    ! Then: All coefficients should be zero
    
    call coeff_rk_gauss(n, a, b, c)
    
    if (any(abs(a) > tol) .or. any(abs(b) > tol) .or. any(abs(c) > tol)) then
      print *, "ERROR: Unsupported Gauss stage count should give zero coefficients"
      errors = errors + 1
    end if
    
  end subroutine test_gauss_unsupported
  
  subroutine test_rk_lobatto_coefficients(errors)
    integer, intent(inout) :: errors
    
    print *, "Testing Runge-Kutta Lobatto coefficients..."
    
    ! Test 3-stage Lobatto method
    call test_lobatto_n3(errors)
    
    ! Test unsupported stage count
    call test_lobatto_unsupported(errors)
    
    if (errors == 0) then
      print *, "  RK Lobatto coefficients test PASSED"
    end if
    
  end subroutine test_rk_lobatto_coefficients
  
  subroutine test_lobatto_n3(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 3
    real(dp) :: a(n,n), ahat(n,n), b(n), c(n)
    real(dp), parameter :: tol = 1.0d-14
    
    ! Given: A 3-stage Lobatto method
    ! When: We compute the coefficients
    ! Then: They should satisfy Lobatto method properties
    
    call coeff_rk_lobatto(n, a, ahat, b, c)
    
    ! Check that b coefficients sum to 1
    if (abs(sum(b) - 1.0d0) > tol) then
      print *, "ERROR: 3-stage Lobatto b coefficients should sum to 1"
      print *, "Sum:", sum(b)
      errors = errors + 1
    end if
    
    ! Check Lobatto property: c(1) = 0, c(n) = 1
    if (abs(c(1)) > tol) then
      print *, "ERROR: Lobatto c(1) should be 0"
      print *, "Got:", c(1)
      errors = errors + 1
    end if
    
    if (abs(c(3) - 1.0d0) > tol) then
      print *, "ERROR: Lobatto c(3) should be 1"
      print *, "Got:", c(3)
      errors = errors + 1
    end if
    
    ! Check that first row of a is zero (Lobatto IIIA property)
    if (any(abs(a(1,:)) > tol)) then
      print *, "ERROR: First row of Lobatto a matrix should be zero"
      errors = errors + 1
    end if
    
    ! Check symmetry of b coefficients for 3-stage
    if (abs(b(1) - b(3)) > tol) then
      print *, "ERROR: Lobatto b(1) should equal b(3)"
      errors = errors + 1
    end if
    
  end subroutine test_lobatto_n3
  
  subroutine test_lobatto_unsupported(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 4  ! Unsupported stage count for Lobatto
    real(dp) :: a(n,n), ahat(n,n), b(n), c(n)
    real(dp), parameter :: tol = 1.0d-14
    
    ! Given: An unsupported stage count for Lobatto
    ! When: We call coeff_rk_lobatto
    ! Then: Arrays remain uninitialized (implementation specific behavior)
    
    ! Initialize arrays to zero before the call
    a = 0.0_dp
    ahat = 0.0_dp
    b = 0.0_dp
    c = 0.0_dp
    
    call coeff_rk_lobatto(n, a, ahat, b, c)
    
    ! Since only n=3 is supported, arrays should remain zero for n=4
    if (any(abs(a) > tol) .or. any(abs(ahat) > tol) .or. &
        any(abs(b) > tol) .or. any(abs(c) > tol)) then
      print *, "ERROR: Unsupported Lobatto stage count should leave coefficients unchanged"
      errors = errors + 1
    end if
    
  end subroutine test_lobatto_unsupported
  
  subroutine test_symplectic_integrator_type(errors)
    integer, intent(inout) :: errors
    type(symplectic_integrator_t) :: si
    type(multistage_integrator_t) :: mi
    
    print *, "Testing symplectic_integrator_t type..."
    
    ! Given: The symplectic_integrator_t and multistage_integrator_t types
    ! When: We initialize them with default values
    ! Then: They should have the expected structure
    
    ! Test symplectic_integrator_t initialization
    si%atol = 1.0d-10
    si%rtol = 1.0d-8
    si%z = [1.0_dp, 0.0_dp, 0.0_dp, 0.1_dp]
    si%pthold = 0.0_dp
    si%ntau = 1000
    si%dt = 1.0d-3
    si%pabs = 0.1_dp
    
    ! Basic checks on data integrity
    if (si%atol /= 1.0d-10) then
      print *, "ERROR: symplectic_integrator_t atol assignment failed"
      errors = errors + 1
    end if
    
    if (size(si%z) /= 4) then
      print *, "ERROR: symplectic_integrator_t z should have 4 components"
      errors = errors + 1
    end if
    
    ! Test multistage_integrator_t initialization
    mi%s = 3
    if (mi%s /= 3) then
      print *, "ERROR: multistage_integrator_t s assignment failed"
      errors = errors + 1
    end if
    
    if (size(mi%alpha) /= S_MAX) then
      print *, "ERROR: multistage_integrator_t alpha array size incorrect"
      errors = errors + 1
    end if
    
    if (size(mi%stages) /= 2*S_MAX) then
      print *, "ERROR: multistage_integrator_t stages array size incorrect"
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  symplectic_integrator_t type test PASSED"
    end if
    
  end subroutine test_symplectic_integrator_type

end program test_orbit_symplectic_base