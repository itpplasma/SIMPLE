!
  subroutine vmecin(rmn,zmn,almn,aiota,phi,sps,axm,axn,s, &
                    nsurfb,nstrb,kparb,flux)
!
  use new_vmec_stuff_mod, only : netcdffile
  use nctools_module
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
  double precision, parameter :: fac_b=1d4,fac_r=1d2
!
  integer :: nsurfb,nstrb,kparb,ncid,i
  double precision :: flux
  double precision, dimension(nstrb)         :: axm,axn
  double precision, dimension(0:kparb)       :: sps,aiota,phi,s
  double precision, dimension(nstrb,0:kparb) :: rmn,zmn,almn
!
  do i=0,kparb
    sps(i)=dble(i)
  enddo
!
  call nc_open(netcdffile, ncid)
!
  call nc_get(ncid, 'phi', phi)
!
  flux=phi(kparb)
  flux=flux*fac_b*fac_r**2
  phi=phi*fac_b*fac_r**2
  s=phi/flux
!
  call nc_get(ncid, 'xm', axm)
  call nc_get(ncid, 'xn', axn)
  call nc_get(ncid, 'iotaf', aiota)
  call nc_get(ncid, 'rmnc', rmn)
  call nc_get(ncid, 'zmns', zmn)
  call nc_get(ncid, 'lmns', almn)
!
  rmn=rmn*fac_r
  zmn=zmn*fac_r
!
  call nc_close(ncid)
!
  end subroutine vmecin
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)
!
  use new_vmec_stuff_mod, only : nper,rmajor
!
  implicit none
! 
  integer :: L1i
  double precision :: RT0,R0i,cbfi,bz0i,bf0
!
  L1i=nper
  RT0=rmajor*1.d2
!
  end subroutine stevvo
