module orbit_symplectic_quasi

use field_can_mod
use orbit_symplectic

implicit none
save

type(FieldCan) :: f
type(SymplecticIntegrator) :: si
!$omp threadprivate(f, si)

contains

subroutine f_midpoint_quasi(n, x, fvec, iflag)
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)  ! = (rend, thend, phend, pphend, rmid)
    double precision, intent(out) :: fvec(n)
    integer, intent(in) :: iflag

    call f_midpoint_part1(si, f, n, x, fvec, iflag)
    call f_midpoint_part2(si, f, n, x, fvec, iflag)
end subroutine f_midpoint_quasi


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine timestep_midpoint_quasi(ierr)
    
      integer, intent(out) :: ierr
    
      integer, parameter :: n = 5
      integer, parameter :: maxit = 256
    
      double precision, dimension(n) :: x
      double precision :: fvec(n)
      integer :: k, ktau, info
    
      ierr = 0
      ktau = 0
      do while(ktau .lt. si%ntau)
        si%pthold = f%pth
    
        x(1:4) = si%z
        x(5) = si%z(1)
    
        call hybrd1(f_midpoint_quasi, n, x, fvec, si%rtol, info)
    
        if (x(1) > 1.0) then
          ierr = 1
          return
        end if
    
        if (x(1) < 0.0) then
          print *, 'r<0, z = ', x(1), si%z(2), si%z(3), x(2)
          x(1) = 0.01
        end if
    
        si%z = x(1:4)
    
        si%kt = si%kt+1
        ktau = ktau+1
      enddo
    
    end subroutine timestep_midpoint_quasi

end module orbit_symplectic_quasi
