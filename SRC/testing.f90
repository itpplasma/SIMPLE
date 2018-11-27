  subroutine testing
!
  use vmec_stuff_mod, only : netcdffile
!
!  implicit none
  implicit real*8 (a-h,o-z),integer(i-n)
!
  logical :: fullset
  integer :: is,npoi
  double precision :: pi
  double precision :: s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,      &
                      R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
  double precision :: ds
  double precision, dimension(3,3) :: cmat,gV,g
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
  double precision, dimension(:), allocatable :: dummy,dummy1
  double precision, dimension(:,:), allocatable :: dummy2d
  double precision, dimension(:,:), allocatable :: dummy2d1
!
  pi=4.d0*atan2(1.d0,1.d0)
!
  npoi=1000
  ds=1.d0/npoi
  theta=0.1d0
  varphi=0.1d0
  vartheta_c=0.5d0
  varphi_c=0.5d0
!
if(.false.) then
!if(.true.) then
  do is=1,npoi
    s=ds*is
!    call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
!                    sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
!                    Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
!    write(1001,*) s,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota
!    write(1001,*) s,alam,dl_dt
!
!    call splint_iota(s,aiota,daiota_ds)
!    call splint_lambda(s,theta,varphi,alam,dl_dt)
!
!    write(1002,*) s,aiota,daiota_ds
!    write(1002,*) s,alam,dl_dt
!
    fullset=.false.
!
    call splint_can_coord(s,vartheta_c,varphi_c,                                           &
                          A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,                 &
                          sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                          B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                          B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                          fullset,G_c)
!
    write (1001,*) s,sqg_c,B_vartheta_c,B_varphi_c
  enddo
stop
endif
!
!
  open(1,file='alpha_lifetime_m.inp')
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*) netcdffile        !name of VMEC file in NETCDF format <=2017 NEW
  close(1)
!
  call allocate_vmec_stuff
!
if(.false.) then
!if(.true.) then
  nt=1000
  np=1000
  allocate(dummy(0:np),dummy1(0:np))
  allocate(dummy2d(12,0:np))
  allocate(dummy2d1(12,0:np))
!  s=3.1415d0/7.d0
  s=0.5d0
  soi=s*98.d0
  do i=0,nt
    teti=2.d0*pi*dfloat(i)/dfloat(nt)
    do j=0,np
      fii=2.d0*pi*dfloat(j)/dfloat(np*5)
!
      vartheta_c=teti
      varphi_c=fii
!
      call can_to_vmec(s,vartheta_c,varphi_c,theta,varphi)
!
!
      call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                      sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                      Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
!print *,s,theta,varphi,dA_theta_ds,dA_phi_ds
!print *,sqg,Bcovar_vartheta*Bctrvr_vartheta+Bcovar_varphi*Bctrvr_varphi
!
      fullset=.false.
!
      call splint_can_coord(s,vartheta_c,varphi_c,                                           &
                            A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,                 &
                            sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                            B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                            B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                            fullset,G_c)
!
!print *,s,vartheta_c,varphi_c,dA_theta_dr,dA_phi_dr
!print *,sqg_c,(-dA_phi_dr*B_vartheta_c+dA_theta_dr*B_varphi_c)/sqg_c
!print *,' '
!      dummy1(j)=Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi
!      dummy1(j)=cjac*dZ_dt*Bctrvr_vartheta+(dZ_dp-dZ_dt*cjac*dl_dp)*Bctrvr_varphi
!
!      call splint_vmec_data(s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,      &
!                            R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
      x(1)=s
      x(2)=theta
      x(3)=varphi
!
      call magfie_vmec(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      dummy2d(1,j)=bmod
      dummy2d(2,j)=hcurl(1)
!
      x(1)=s
      x(2)=vartheta_c
      x(3)=varphi_c
!
      call magfie_can(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      dummy2d1(1,j)=bmod
      dummy2d1(2,j)=hcurl(1)
    enddo
    do j=1,2
      write(2000+j,*) dummy2d(j,:)
      write(2100+j,*) dummy2d1(j,:)
    enddo
  enddo
stop
!
!  call deallocate_vmec_stuff
endif
!
  end subroutine testing
