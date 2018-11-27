!
!  module chamb_mod
!    double precision :: rbig,rcham2
!  end module chamb_mod
!
!  module parmot_mod
!    double precision :: rmu,ro0,eeff
!  end module parmot_mod
!
  module vmec_stuff_mod
    character*32     :: netcdffile
    integer          :: nsurfm,nstrm,nper,kpar
    double precision :: rmajor
    double precision, dimension(:),     allocatable :: axm,axn,soa
    double precision, dimension(:,:),   allocatable :: saiota,ss,ssps,sphi
    double precision, dimension(:,:,:), allocatable :: slmn,srmn,szmn
  end module vmec_stuff_mod
!
!module gbpi_mod
!  character*24 :: filed
!  integer :: ierrfield
!end module gbpi_mod
