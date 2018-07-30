  subroutine testing
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
  theta=0.5d0
  varphi=0.5d0
  vartheta_c=0.5d0
  varphi_c=0.5d0
!
!if(.false.) then
if(.true.) then
  do is=1,npoi
    s=ds*is
    call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                    sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                    Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
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
    write (1001,*) s,sqg,sqg_c,Bcovar_vartheta,B_vartheta_c,Bcovar_varphi,B_varphi_c
    write (1002,*) s,sqg,Bcovar_vartheta,Bcovar_varphi,alam,dl_ds
  enddo
stop
endif
!
end
