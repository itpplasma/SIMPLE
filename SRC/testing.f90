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
  vartheta_c=0.5d0
  varphi_c=0.5d0
!
!if(.false.) then
if(.true.) then
  do is=1,npoi
    s=ds*is
    r=s
    call splint_can_coord_moreders(r,vartheta_c,varphi_c,                                           &
                                   A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,                 &
                                   sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                                   B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                                   B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                                   d2sqg_c_dr2,d2B_vartheta_c_dr2,d2B_varphi_c_dr2,d3A_phi_dr3)
!
    Bctr_vartheta=-dA_phi_dr/sqg_c
    Bctr_varphi=dA_theta_dr/sqg_c
!
    bmod2=Bctr_vartheta*B_vartheta_c+Bctr_varphi*B_varphi_c
    bmod=sqrt(bmod2)
!
    dbmod_dr=0.5d0*((dA_theta_dr*dB_varphi_c_dr-dA_phi_dr*dB_vartheta_c_dr-d2A_phi_dr2*B_vartheta_c) &
             /bmod2-dsqg_c_dr)*bmod/sqg_c
!
    d2bmod_dr2=(dA_theta_dr*d2B_varphi_c_dr2-dA_phi_dr*d2B_vartheta_c_dr2                            &
              - 2.d0*d2A_phi_dr2*dB_vartheta_c_dr-d3A_phi_dr3*B_vartheta_c                           &
              - bmod2*d2sqg_c_dr2)/(2.d0*bmod*sqg_c)-dbmod_dr**2/bmod-2.d0*dbmod_dr*dsqg_c_dr/sqg_c
!
    h_theta=B_vartheta_c/bmod
    dh_theta_dr=(dB_vartheta_c_dr-h_theta*dbmod_dr)/bmod
    d2h_theta_dr2=(d2B_vartheta_c_dr2-2.d0*dh_theta_dr*dbmod_dr-h_theta*d2bmod_dr2)/bmod
!
    h_phi=B_varphi_c/bmod
    dh_phi_dr=(dB_varphi_c_dr-h_phi*dbmod_dr)/bmod
    d2h_phi_dr2=(d2B_varphi_c_dr2-2.d0*dh_phi_dr*dbmod_dr-h_phi*d2bmod_dr2)/bmod
!
    write (1001,*) s,dbmod_dr,d2bmod_dr2
    write (1002,*) s,dh_theta_dr,d2h_theta_dr2
    write (1003,*) s,dh_phi_dr,d2h_phi_dr2
  enddo
stop
endif
!
end
