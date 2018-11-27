!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine chamb_can(y,phi,ierr)
!
! checks whether the point is inside the vacuum chamber
!  Input parameters:
!            formal: y(i) - coordinates on the poloidal cross section
!                    phi  - toroidal angle
! Outout parameters:
!            formal: ierr -error code (0 if the point is inside 1 - othervice)
!
      implicit none
!
      integer :: ierr
      double precision :: phi
      double precision, dimension(2) :: y
!
      if(y(1).ge.1.d0) then
        ierr = 1
      else
        ierr = 0
      endif
!
      end subroutine chamb_can
