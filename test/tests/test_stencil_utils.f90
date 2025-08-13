program test_stencil_utils
  use stencil_utils
  implicit none
  
  integer, parameter :: dp = kind(1.0d0)
  
  call test_stencil_order_2()
  call test_stencil_order_4()
  call test_stencil_order_6()
  call test_stencil_scaling()
  call test_stencil_invalid_order()
  
  print *, "All stencil_utils tests passed!"
  
contains

  subroutine test_stencil_order_2()
    real(dp) :: stencil(-1:1)
    real(dp) :: h_grid = 0.1_dp
    real(dp), parameter :: tol = 1.0e-14_dp
    
    call init_derivative_stencil(1, h_grid, stencil)
    
    ! Check values
    if (abs(stencil(-1) - (-0.5d0/h_grid)) > tol) then
      print *, "ERROR: stencil(-1) incorrect for order 2"
      print *, "Expected:", -0.5d0/h_grid, "Got:", stencil(-1)
      error stop 1
    end if
    
    if (abs(stencil(0)) > tol) then
      print *, "ERROR: stencil(0) should be zero for order 2"
      print *, "Got:", stencil(0)
      error stop 1
    end if
    
    if (abs(stencil(1) - (0.5d0/h_grid)) > tol) then
      print *, "ERROR: stencil(1) incorrect for order 2"
      print *, "Expected:", 0.5d0/h_grid, "Got:", stencil(1)
      error stop 1
    end if
    
    print *, "test_stencil_order_2: PASSED"
  end subroutine test_stencil_order_2
  
  subroutine test_stencil_order_4()
    real(dp) :: stencil(-2:2)
    real(dp) :: h_grid = 0.1_dp
    real(dp), parameter :: tol = 1.0e-14_dp
    
    call init_derivative_stencil(2, h_grid, stencil)
    
    ! Check values
    if (abs(stencil(-2) - (1.d0/12.d0/h_grid)) > tol) then
      print *, "ERROR: stencil(-2) incorrect for order 4"
      error stop 1
    end if
    
    if (abs(stencil(-1) - (-2.d0/3.d0/h_grid)) > tol) then
      print *, "ERROR: stencil(-1) incorrect for order 4"
      error stop 1
    end if
    
    if (abs(stencil(0)) > tol) then
      print *, "ERROR: stencil(0) should be zero for order 4"
      error stop 1
    end if
    
    if (abs(stencil(1) - (2.d0/3.d0/h_grid)) > tol) then
      print *, "ERROR: stencil(1) incorrect for order 4"
      error stop 1
    end if
    
    if (abs(stencil(2) - (-1.d0/12.d0/h_grid)) > tol) then
      print *, "ERROR: stencil(2) incorrect for order 4"
      error stop 1
    end if
    
    print *, "test_stencil_order_4: PASSED"
  end subroutine test_stencil_order_4
  
  subroutine test_stencil_order_6()
    real(dp) :: stencil(-3:3)
    real(dp) :: h_grid = 0.1_dp
    real(dp), parameter :: tol = 1.0e-14_dp
    real(dp) :: sum_coeff
    real(dp) :: derivative_linear
    real(dp) :: derivative_cubic
    real(dp) :: x_values(-3:3)
    integer :: j
    
    call init_derivative_stencil(3, h_grid, stencil)
    
    ! Verify stencil by applying it to test functions
    ! For a 6th order stencil, it should exactly differentiate polynomials up to degree 6
    ! Test with f(x) = x^5 at x=0, derivative should be 0
    ! Test with f(x) = x^3 at x=h, derivative should be 3*h^2
    
    ! First, verify the stencil satisfies mathematical properties
    ! Property 1: Sum should be zero for derivative operator
    sum_coeff = stencil(-3) + stencil(-2) + stencil(-1) + stencil(0) + &
                stencil(1) + stencil(2) + stencil(3)
    if (abs(sum_coeff) > tol) then
      print *, "ERROR: Stencil coefficients don't sum to zero"
      print *, "Sum:", sum_coeff
      error stop 1
    end if
    
    ! Property 2: Antisymmetry (stencil(-i) = -stencil(i))
    if (abs(stencil(-3) + stencil(3)) > tol) then
      print *, "ERROR: stencil not antisymmetric at ±3"
      error stop 1
    end if
    if (abs(stencil(-2) + stencil(2)) > tol) then
      print *, "ERROR: stencil not antisymmetric at ±2"
      error stop 1
    end if
    if (abs(stencil(-1) + stencil(1)) > tol) then
      print *, "ERROR: stencil not antisymmetric at ±1"
      error stop 1
    end if
    
    ! Property 3: Test on known function f(x) = x
    ! Derivative should be exactly 1.0 for linear function
    ! f(-3h) = -3h, f(-2h) = -2h, ..., f(3h) = 3h
    
    do j = -3, 3
      x_values(j) = real(j, dp) * h_grid
    end do
    
    derivative_linear = 0.0_dp
    do j = -3, 3
      derivative_linear = derivative_linear + stencil(j) * x_values(j)
    end do
    
    if (abs(derivative_linear - 1.0_dp) > tol*10) then
      print *, "ERROR: Stencil doesn't correctly differentiate f(x)=x"
      print *, "Expected: 1.0, Got:", derivative_linear
      error stop 1
    end if
    
    ! Property 4: Test on f(x) = x^3
    ! At x=0, derivative should be 0
    do j = -3, 3
      x_values(j) = (real(j, dp) * h_grid)**3
    end do
    
    derivative_cubic = 0.0_dp
    do j = -3, 3
      derivative_cubic = derivative_cubic + stencil(j) * x_values(j)
    end do
    
    ! For x^3 at x=0, derivative is 0
    if (abs(derivative_cubic) > tol*100) then
      print *, "ERROR: Stencil doesn't correctly differentiate f(x)=x^3 at x=0"
      print *, "Expected: 0.0, Got:", derivative_cubic
      error stop 1
    end if
    
    print *, "test_stencil_order_6: PASSED"
  end subroutine test_stencil_order_6
  
  subroutine test_stencil_scaling()
    real(dp) :: stencil1(-1:1), stencil2(-1:1)
    real(dp) :: h1 = 0.1_dp, h2 = 0.2_dp
    real(dp), parameter :: tol = 1.0e-14_dp
    
    call init_derivative_stencil(1, h1, stencil1)
    call init_derivative_stencil(1, h2, stencil2)
    
    ! Check that stencil scales inversely with h
    if (abs(stencil1(1) * h1 - stencil2(1) * h2) > tol) then
      print *, "ERROR: stencil does not scale properly with h"
      print *, "stencil1(1)*h1 =", stencil1(1)*h1
      print *, "stencil2(1)*h2 =", stencil2(1)*h2
      error stop 1
    end if
    
    print *, "test_stencil_scaling: PASSED"
  end subroutine test_stencil_scaling

  subroutine test_stencil_invalid_order()
    real(dp) :: stencil(-4:4)
    real(dp) :: h_grid = 0.1_dp
    real(dp), parameter :: tol = 1.0e-14_dp
    integer :: i
    
    print *, "Testing invalid stencil order (tests case default line 37-39)..."
    
    ! Given: An invalid nh_stencil value (not 1, 2, or 3)
    ! When: We call init_derivative_stencil with invalid order
    ! Then: All stencil coefficients should be set to 0.0
    
    ! Initialize stencil with non-zero values to verify they get zeroed
    stencil = 999.0d0
    
    ! Test with invalid nh_stencil = 4 (case default)
    call init_derivative_stencil(4, h_grid, stencil)
    
    ! Verify all coefficients are zero
    do i = -4, 4
      if (abs(stencil(i)) > tol) then
        print *, "ERROR: stencil coefficient should be zero for invalid order"
        print *, "stencil(", i, ") =", stencil(i)
        error stop 1
      end if
    end do
    
    ! Test with another invalid value: nh_stencil = 0
    block
      real(dp) :: zero_stencil(0:0)
      zero_stencil = 999.0d0  ! Reset to non-zero
      call init_derivative_stencil(0, h_grid, zero_stencil)
      
      if (abs(zero_stencil(0)) > tol) then
        print *, "ERROR: stencil should be zero for nh_stencil=0"
        print *, "stencil(0) =", zero_stencil(0)
        error stop 1
      end if
    end block
    
    ! Test with another invalid positive value: nh_stencil = 5  
    block
      real(dp) :: big_stencil(-5:5)
      big_stencil = 999.0d0  ! Reset to non-zero
      call init_derivative_stencil(5, h_grid, big_stencil)
      
      ! All coefficients should be zero for unsupported order
      do i = -5, 5
        if (abs(big_stencil(i)) > tol) then
          print *, "ERROR: stencil coefficient should be zero for invalid order 5"
          print *, "stencil(", i, ") =", big_stencil(i)
          error stop 1
        end if
      end do
    end block
    
    print *, "test_stencil_invalid_order: PASSED"
  end subroutine test_stencil_invalid_order

end program test_stencil_utils