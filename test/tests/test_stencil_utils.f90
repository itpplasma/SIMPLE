program test_stencil_utils
  use stencil_utils
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  call test_stencil_polynomial_exactness()
  call test_stencil_invalid_order()

  print *, "All stencil_utils tests passed!"

contains

  subroutine test_stencil_polynomial_exactness()
    ! A central stencil of half-width nh has order 2*nh: applied to any
    ! polynomial of degree <= 2*nh it reproduces the exact first derivative.
    real(dp), parameter :: h_values(2) = [0.1_dp, 0.037_dp]
    real(dp), parameter :: x0 = 0.3_dp
    real(dp), parameter :: tol = 1.0e-12_dp
    real(dp) :: stencil(-3:3), deriv, exact
    integer :: nh, ih, k, j

    do nh = 1, 3
      do ih = 1, size(h_values)
        stencil = 0.0_dp
        call init_derivative_stencil(nh, h_values(ih), stencil(-nh:nh))
        do k = 0, 2*nh
          deriv = 0.0_dp
          do j = -nh, nh
            deriv = deriv + stencil(j)*(x0 + real(j, dp)*h_values(ih))**k
          end do
          if (k == 0) then
            exact = 0.0_dp
          else
            exact = real(k, dp)*x0**(k - 1)
          end if
          if (abs(deriv - exact) > tol) then
            print *, "ERROR: stencil not exact for polynomial degree", k
            print *, "nh =", nh, "h =", h_values(ih)
            print *, "Expected:", exact, "Got:", deriv
            error stop 1
          end if
        end do
      end do
    end do

    print *, "test_stencil_polynomial_exactness: PASSED"
  end subroutine test_stencil_polynomial_exactness

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