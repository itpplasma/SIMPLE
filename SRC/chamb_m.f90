!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chamb(y,phi,ierr)
!
! checks whether the point is inside the vacuum chamber
!  Input parameters:
!            formal: y(i) - coordinates on the poloidal cross section
!                    phi  - toroidal angle
! Outout parameters:
!            formal: ierr -error code (0 if the point is inside 1 - othervice)
!
      use chamb_mod
      use vmec_stuff_mod, only : kpar                                    !<=2017 NEW
!
      implicit none
!
!<=2017 NEW      integer :: nsurfm,nstrm,nper,kpar !<=2017
!<=2017 NEW      include 'parvmec.f'               !<=2017
      integer :: ierr
      double precision y,phi
!
      dimension y(2)
!
! 27.06.2010
      logical first
      integer :: nier
      save nier,aldd22,first
! 02.07.2010      double precision rbig0
! 23.09.2013      double precision rbig0,facr
      double precision facr,aldd2,aldd22
      data first/.true./
      if(first) then
! 11.11.2012      facr=0.1d0
      facr=1.0d0
      open(92,form='FORMATTED',file='chamb0.dat')
      nier=0
! 11.11.2012      rbig0=93.0/facr
! 23.09.2013      rbig0=100.0d0/facr
      aldd2=314.159265/facr
      aldd22=aldd2**2
      rcham2=(kpar-1)**2 !<=2017
print *,rcham2,kpar-1
      first=.false.
      endif

! 27.06.201 end
! 27.06.2010      if((y(1)-rbig)**2+y(2)**2.lt.rcham2) then
! 23.09.2013      if((y(1)-rbig0)**2+y(2)**2.lt.rcham2) then
! 23.09.2013        ierr=0
! 23.09.2013      else
! 23.09.2013        ierr=1
! 23.09.2013      endif
! 23.09.2013

! 05.03.2014      if(y(1)**2.ge.rcham2.or.y(2)**2.ge.aldd22) then
! 08.03.2014      if(y(1).ge.70.) then
!      if(y(1)**2+y(2)**2.ge.4900.) then
!      if(y(1)**2+y(2)**2.ge.8100.) then
! 07.04.2016      if(y(1)**2+y(2)**2.ge.10000.) then
! 13.07.2016      if(y(1)**2+y(2)**2.ge.5329.) then
!2017=>      if(y(1)**2+y(2)**2.ge.9409.) then
      if(y(1)**2+y(2)**2.ge.rcham2) then  !<=2017
        ierr=1
      else
        ierr=0
      endif
! 23.09.2013 end
      if(ierr.ne.0) then
      nier=nier+1
      print*,y(1),y(2),nier,'=nier, phi=',phi
      write(92,*)y(1),y(2),nier,'=nier, phi=',phi
      endif
      return
      end
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
