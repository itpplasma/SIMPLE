module array_utils
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  
  private
  public :: init_derivative_factors
  
contains

  !> Initialize factorial-based derivative factors for polynomial derivatives
  !> derf1(k) = (k-1)
  !> derf2(k) = (k-1)*(k-2) 
  !> derf3(k) = (k-1)*(k-2)*(k-3)
  pure subroutine init_derivative_factors(ns_max, derf1, derf2, derf3)
    integer, intent(in) :: ns_max
    double precision, intent(out) :: derf1(ns_max), derf2(ns_max), derf3(ns_max)
    integer :: k
    
    do k = 1, ns_max
      derf1(k) = dble(k-1)
      derf2(k) = dble((k-1)*(k-2))
      derf3(k) = dble((k-1)*(k-2)*(k-3))
    enddo
    
  end subroutine init_derivative_factors

end module array_utils