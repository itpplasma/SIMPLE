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
  npoi=1000
  ds=1.d0/npoi
  vartheta_c=0.5d0
  varphi_c=0.5d0
!
print *,'start test'
do ii=1,100000
if(ii/1000*1000.eq.ii) print *,ii
  do is=1,npoi
    s=ds*is
    r=s
    mode_secders=2
!
    call splint_can_coord(fullset,mode_secders,r,vartheta_c,varphi_c,                      &
                          A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3,     &
                          sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                          B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                          B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                          d2sqg_rr,d2sqg_rt,d2sqg_rp,d2sqg_tt,d2sqg_tp,d2sqg_pp,           &
                          d2bth_rr,d2bth_rt,d2bth_rp,d2bth_tt,d2bth_tp,d2bth_pp,           &
                          d2bph_rr,d2bph_rt,d2bph_rp,d2bph_tt,d2bph_tp,d2bph_pp,G_c)
!
    bmod2=(B_varphi_c*dA_theta_dr-B_vartheta_c*dA_phi_dr)/sqg_c
    dbmod2_dr=(dB_varphi_c_dr*dA_theta_dr-dB_vartheta_c_dr*dA_phi_dr-B_vartheta_c*d2A_phi_dr2-bmod2*dsqg_c_dr)/sqg_c
    dbmod2_dt=(dB_varphi_c_dt*dA_theta_dr-dB_vartheta_c_dt*dA_phi_dr-bmod2*dsqg_c_dt)/sqg_c
    dbmod2_dp=(dB_varphi_c_dp*dA_theta_dr-dB_vartheta_c_dp*dA_phi_dr-bmod2*dsqg_c_dp)/sqg_c
!
    bmod=sqrt(bmod2)
    twobmod=2.d0*bmod
    dbmod_dr=dbmod2_dr/twobmod
    dbmod_dt=dbmod2_dt/twobmod
    dbmod_dp=dbmod2_dp/twobmod
!
    if(mode_secders.eq.1) then
      d2bmod2_rr=(d2bph_rr*dA_theta_dr-d2bth_rr*dA_phi_dr-2.d0*dB_vartheta_c_dr*d2A_phi_dr2-B_vartheta_c*d3A_phi_dr3 &
                 - 2.d0*dsqg_c_dr*dbmod2_dr-bmod2*d2sqg_rr)/sqg_c 
!
      d2bmod_rr=d2bmod2_rr/twobmod-dbmod_dr**2/bmod
    elseif(mode_secders.eq.2) then
      d2bmod2_rr=(d2bph_rr*dA_theta_dr-d2bth_rr*dA_phi_dr-2.d0*dB_vartheta_c_dr*d2A_phi_dr2-B_vartheta_c*d3A_phi_dr3 &
                 - 2.d0*dsqg_c_dr*dbmod2_dr-bmod2*d2sqg_rr)/sqg_c 
      d2bmod2_rt=(d2bph_rt*dA_theta_dr-d2bth_rt*dA_phi_dr-dB_vartheta_c_dt*d2A_phi_dr2                               &
                 - dsqg_c_dr*dbmod2_dt-dsqg_c_dt*dbmod2_dr-bmod2*d2sqg_rt)/sqg_c 
      d2bmod2_rp=(d2bph_rp*dA_theta_dr-d2bth_rp*dA_phi_dr-dB_vartheta_c_dp*d2A_phi_dr2                               &
                 - dsqg_c_dr*dbmod2_dp-dsqg_c_dp*dbmod2_dr-bmod2*d2sqg_rp)/sqg_c 
      d2bmod2_tt=(d2bph_tt*dA_theta_dr-d2bth_tt*dA_phi_dr-2.d0*dsqg_c_dt*dbmod2_dt-bmod2*d2sqg_tt)/sqg_c 
      d2bmod2_tp=(d2bph_tp*dA_theta_dr-d2bth_tp*dA_phi_dr-dsqg_c_dt*dbmod2_dp-dsqg_c_dp*dbmod2_dt-bmod2*d2sqg_tp)/sqg_c 
      d2bmod2_pp=(d2bph_pp*dA_theta_dr-d2bth_pp*dA_phi_dr-2.d0*dsqg_c_dp*dbmod2_dp-bmod2*d2sqg_pp)/sqg_c 
!
      d2bmod_rr=d2bmod2_rr/twobmod-dbmod_dr**2/bmod
      d2bmod_rt=d2bmod2_rt/twobmod-dbmod_dr*dbmod_dt/bmod
      d2bmod_rp=d2bmod2_rp/twobmod-dbmod_dr*dbmod_dp/bmod
      d2bmod_tt=d2bmod2_tt/twobmod-dbmod_dt**2/bmod
      d2bmod_tp=d2bmod2_tp/twobmod-dbmod_dt*dbmod_dp/bmod
      d2bmod_pp=d2bmod2_pp/twobmod-dbmod_dp**2/bmod
    endif
!
  enddo
enddo
stop
!
end
