  subroutine testing
!
!  implicit none
  implicit real*8 (a-h,o-z),integer(i-n)
!
  logical :: fullset=.false.
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
  npoi=100000
  ds=1.d0/npoi
!npoi=100
  vartheta_c=0.5d0
  varphi_c=0.5d0
!
ntheta=100
allocate(dummy(ntheta))
do is=1,npoi
!  do it=1,ntheta
    s=ds*is
    r=s
    mode_secders=2
!    vartheta_c=2.d0*pi*dble(it)/dble(ntheta)
  call vmec_field(s,vartheta_c,varphi_c,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
write(103,*) s,Bcovar_r,Bcovar_vartheta,Bcovar_varphi,sqg
!
if(.false.) then
    call splint_can_coord(fullset,mode_secders,r,vartheta_c,varphi_c,                      &
                          A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3,     &
                          sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                          B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                          B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                          d2sqg_rr,d2sqg_rt,d2sqg_rp,d2sqg_tt,d2sqg_tp,d2sqg_pp,           &
                          d2bth_rr,d2bth_rt,d2bth_rp,d2bth_tt,d2bth_tp,d2bth_pp,           &
                          d2bph_rr,d2bph_rt,d2bph_rp,d2bph_tt,d2bph_tp,d2bph_pp,G_c)
!
!    dummy(it)=B_vartheta_c
!  enddo
!  write(101,*) dummy
write(102,*) r,B_vartheta_c,B_varphi_c,sqg_c
endif
enddo
stop
!
end
