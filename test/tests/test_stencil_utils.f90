program test_stencil_utils
  use stencil_utils
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  
  call test_stencil_order_2()
  call test_stencil_order_4()
  call test_stencil_order_6()
  call test_stencil_scaling()
  call test_stencil_invalid_order()
  
  print *, "All stencil_utils tests passed!"
  
contains

  subroutine test_stencil_order_2()
    double precision :: stencil(-1:1)
    double precision :: h_grid = 0.1d0
    double precision, parameter :: tol = 1.0d-14
    
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
    double precision :: stencil(-2:2)
    double precision :: h_grid = 0.1d0
    double precision, parameter :: tol = 1.0d-14
    
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
    double precision :: stencil(-3:3)
    double precision :: h_grid = 0.1d0
    double precision, parameter :: tol = 1.0d-14
    
    call init_derivative_stencil(3, h_grid, stencil)
    
    ! Check symmetry
    if (abs(stencil(-3) + stencil(3)) > tol) then
      print *, "ERROR: stencil not antisymmetric for order 6"
      error stop 1
    end if
    
    if (abs(stencil(-2) + stencil(2)) > tol) then
      print *, "ERROR: stencil not antisymmetric for order 6"
      error stop 1
    end if
    
    if (abs(stencil(-1) + stencil(1)) > tol) then
      print *, "ERROR: stencil not antisymmetric for order 6"
      error stop 1
    end if
    
    print *, "test_stencil_order_6: PASSED"
  end subroutine test_stencil_order_6
  
  subroutine test_stencil_scaling()
    double precision :: stencil1(-1:1), stencil2(-1:1)
    double precision :: h1 = 0.1d0, h2 = 0.2d0
    double precision, parameter :: tol = 1.0d-14
    
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
    double precision :: stencil(-4:4)
    double precision :: h_grid = 0.1d0
    double precision, parameter :: tol = 1.0d-14
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