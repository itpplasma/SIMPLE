module stencil_utils
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  
  private
  public :: init_derivative_stencil
  
contains

  pure subroutine init_derivative_stencil(nh_stencil, h_grid, stencil)
    integer, intent(in) :: nh_stencil
    double precision, intent(in) :: h_grid
    double precision, intent(out) :: stencil(-nh_stencil:nh_stencil)
    
    select case(nh_stencil)
    case(1)
      ! 2nd order centered difference
      stencil(-1) = -0.5d0
      stencil(0) = 0.0d0
      stencil(1) = 0.5d0
    case(2)
      ! 4th order centered difference
      stencil(-2) = 1.d0/12.d0
      stencil(-1) = -2.d0/3.d0
      stencil(0) = 0.0d0
      stencil(1) = 2.d0/3.d0
      stencil(2) = -1.d0/12.d0
    case(3)
      ! 6th order centered difference
      stencil(-3) = -1.d0/60.d0
      stencil(-2) = 0.15d0
      stencil(-1) = -0.75d0
      stencil(0) = 0.0d0
      stencil(1) = 0.75d0
      stencil(2) = -0.15d0
      stencil(3) = 1.d0/60.d0
    case default
      ! This should never happen if input validation is done properly
      stencil = 0.0d0
    end select
    
    ! Scale by grid spacing
    stencil = stencil / h_grid
    
  end subroutine init_derivative_stencil

end module stencil_utils