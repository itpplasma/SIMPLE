!
  subroutine new_allocate_vmec_stuff
!
  use new_vmec_stuff_mod
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
  rmajor=rmajor*vmec_RZ_scale
!
  call nc_close(ncid)
!
  allocate(axm(nstrm),axn(nstrm),soa(0:kpar))
  allocate(aiota(0:kpar),sps(0:kpar),phi(0:kpar),s(0:kpar))
  allocate(rmnc(nstrm,0:kpar),zmnc(nstrm,0:kpar),almnc(nstrm,0:kpar))
  allocate(rmns(nstrm,0:kpar),zmns(nstrm,0:kpar),almns(nstrm,0:kpar))
!
  end subroutine new_allocate_vmec_stuff
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine new_deallocate_vmec_stuff
!
  use new_vmec_stuff_mod
!
  implicit none
!
  deallocate(axm,axn,soa,aiota,sps,phi,s,rmnc,zmnc,almnc,rmns,zmns,almns)
!
  end subroutine new_deallocate_vmec_stuff
