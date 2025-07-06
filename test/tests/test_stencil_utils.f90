program test_stencil_utils
  use stencil_utils
  implicit none
  
  call test_stencil_order_2()
  call test_stencil_order_4()
  call test_stencil_order_6()
  call test_stencil_scaling()
  
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

end program test_stencil_utils