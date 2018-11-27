!
  subroutine allocate_vmec_stuff
!
  use vmec_stuff_mod
  use nctools_module
! 
  implicit none
!
  integer :: ncid
  integer, dimension(2) :: lens
!
  call nc_open(netcdffile, ncid)
!
  call nc_inq_dim(ncid, 'lmns', lens)
!
  call nc_get(ncid, 'nfp', nper)
!
  call nc_get(ncid, 'Rmajor_p', rmajor)
!
  nstrm=lens(1)
  nsurfm=lens(2)
  kpar=nsurfm-1
!
  call nc_close(ncid)
!
  allocate(axm(nstrm),axn(nstrm),soa(0:kpar))
  allocate(saiota(4,kpar),ssps(4,kpar),sphi(4,kpar),ss(4,kpar))
  allocate(srmn(4,nstrm,kpar),szmn(4,nstrm,kpar),slmn(4,nstrm,kpar))
!
  end subroutine allocate_vmec_stuff
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine deallocate_vmec_stuff
!
  use vmec_stuff_mod
!
  implicit none
!
  deallocate(axm,axn,soa,saiota,ssps,sphi,ss,srmn,szmn,slmn)
!
  end subroutine deallocate_vmec_stuff
