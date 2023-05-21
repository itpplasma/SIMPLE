!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine magfie_vmec(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
! Computes magnetic field module in units of the magnetic code  - bmod,
! square root of determinant of the metric tensor               - sqrtg,
! derivatives of the logarythm of the magnetic field module
! over coordinates                                              - bder,
! covariant componets of the unit vector of the magnetic
! field direction                                               - hcovar,
! contravariant components of this vector                       - hctrvr,
! contravariant component of the curl of this vector            - hcurl
! Order of coordinates is the following: x(1)=s (normalized toroidal flux),
! x(2)=theta (VMEC poloidal angle), x(3)=varphi (geometrical toroidal angle).
!
!  Input parameters:
!            formal:  x(3)             - array of VMEC coordinates
!  Output parameters:
!            formal:  bmod
!                     sqrtg
!                     bder(3)          - derivatives of $\log(B)$
!                     hcovar(3)        - covariant components of unit vector $\bh$ along $\bB$
!                     hctrvr(3)        - contra-variant components of unit vector $\bh$ along $\bB$
!                     hcurl(3)         - contra-variant components of curl of $\bh$
!
!  Called routines: vmec_field
!
  implicit none
!
  double precision, parameter :: twopi=2.d0*3.14159265358979d0, hs=1.d-3, ht=hs*twopi, hp=ht/5.d0
!
  double precision, intent(out) :: bmod,sqrtg
  double precision :: s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                      sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                      Bcovar_r,Bcovar_vartheta,Bcovar_varphi
  double precision :: cjac,bcov_s_vmec,bcov_t_vmec,bcov_p_vmec
  double precision :: dhs_dt,dhs_dp,dht_ds,dht_dp,dhp_ds,dhp_dt
  double precision, dimension(3), intent(in) :: x
  double precision, dimension(3), intent(out) :: bder,hcovar,hctrvr,hcurl
!
! Begin derivatives over s
!
  theta=x(2)
  varphi=x(3)
  s=x(1)+hs
!
  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(1)=bmod
  dht_ds=bcov_t_vmec/bmod
  dhp_ds=bcov_p_vmec/bmod
!
  s=x(1)-hs
!
  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(1)=(bder(1)-bmod)/(2.d0*hs)
  dht_ds=(dht_ds-bcov_t_vmec/bmod)/(2.d0*hs)
  dhp_ds=(dhp_ds-bcov_p_vmec/bmod)/(2.d0*hs)
!
! End derivatives over s
!
!-------------------------
!
! Begin derivatives over theta
!
  s=x(1)
  theta=x(2)+ht
!
  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(2)=bmod
  dhs_dt=bcov_s_vmec/bmod
  dhp_dt=bcov_p_vmec/bmod
!
  theta=x(2)-ht
!
  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(2)=(bder(2)-bmod)/(2.d0*ht)
  dhs_dt=(dhs_dt-bcov_s_vmec/bmod)/(2.d0*ht)
  dhp_dt=(dhp_dt-bcov_p_vmec/bmod)/(2.d0*ht)
!
! End derivatives over theta
!
!-------------------------
!
! Begin derivatives over varphi
!
  theta=x(2)
  varphi=x(3)+hp
!
  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(3)=bmod
  dhs_dp=bcov_s_vmec/bmod
  dht_dp=bcov_t_vmec/bmod
!
  varphi=x(3)-hp
!
  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  bder(3)=(bder(3)-bmod)/(2.d0*hp)
  dhs_dp=(dhs_dp-bcov_s_vmec/bmod)/(2.d0*hp)
  dht_dp=(dht_dp-bcov_t_vmec/bmod)/(2.d0*hp)
!
! End derivatives over varphi
!
!-------------------------
!
  varphi=x(3)
!
  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
  bmod=sqrt(Bctrvr_vartheta*Bcovar_vartheta+Bctrvr_varphi*Bcovar_varphi)
  cjac=1.d0+dl_dt
  sqrtg=sqg*cjac
  bder=bder/bmod
  bcov_s_vmec=Bcovar_r+dl_ds*Bcovar_vartheta
  bcov_t_vmec=(1.d0+dl_dt)*Bcovar_vartheta
  bcov_p_vmec=Bcovar_varphi+dl_dp*Bcovar_vartheta
  hcovar(1)=bcov_s_vmec/bmod
  hcovar(2)=bcov_t_vmec/bmod
  hcovar(3)=bcov_p_vmec/bmod
  hctrvr(1)=0.d0
  hctrvr(2)=(Bctrvr_vartheta-dl_dp*Bctrvr_varphi)/(cjac*bmod)
  hctrvr(3)=Bctrvr_varphi/bmod
  hcurl(1)=(dhp_dt-dht_dp)/sqrtg
  hcurl(2)=(dhs_dp-dhp_ds)/sqrtg
  hcurl(3)=(dht_ds-dhs_dt)/sqrtg
!
  end subroutine magfie_vmec
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine magfie_can(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
!
! Computes magnetic field module in units of the magnetic code  - bmod,
! square root of determinant of the metric tensor               - sqrtg,
! derivatives of the logarythm of the magnetic field module
! over coordinates                                              - bder,
! covariant componets of the unit vector of the magnetic
! field direction                                               - hcovar,
! contravariant components of this vector                       - hctrvr,
! contravariant component of the curl of this vector            - hcurl
! Order of coordinates is the following: x(1)=s (normalized toroidal flux),
! x(2)=vartheta_c (canonical poloidal angle), x(3)=varphi_c (canonical toroidal angle).
!
!  Input parameters:
!            formal:  x(3)             - array of canonical coordinates
!  Output parameters:
!            formal:  bmod
!                     sqrtg
!                     bder(3)          - derivatives of $\log(B)$
!                     hcovar(3)        - covariant components of unit vector $\bh$ along $\bB$
!                     hctrvr(3)        - contra-variant components of unit vector $\bh$ along $\bB$
!                     hcurl(3)         - contra-variant components of curl of $\bh$
!
!  Called routines: canonical_field
!
  implicit none
!
  logical :: fullset
  integer :: mode_secders
  double precision :: bmod,sqrtg
  double precision :: r,vartheta_c,varphi_c,                                           &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,d3A_phi_dr3,     &
                      sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                      B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                      B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c,     &
                      d2sqg_rr,d2sqg_rt,d2sqg_rp,d2sqg_tt,d2sqg_tp,d2sqg_pp,           &
                      d2bth_rr,d2bth_rt,d2bth_rp,d2bth_tt,d2bth_tp,d2bth_pp,           &
                      d2bph_rr,d2bph_rt,d2bph_rp,d2bph_tt,d2bph_tp,d2bph_pp
  double precision :: Bctr_vartheta,Bctr_varphi,bmod2
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
!
  r=x(1)
  vartheta_c=x(2)
  varphi_c=x(3)
!
  fullset=.false.
  mode_secders=0
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
  sqrtg=sqg_c
!
  Bctr_vartheta=-dA_phi_dr/sqg_c
  Bctr_varphi=dA_theta_dr/sqg_c
!
  bmod2=Bctr_vartheta*B_vartheta_c+Bctr_varphi*B_varphi_c
  bmod=sqrt(bmod2)
!
  bder(1)=0.5d0*((dA_theta_dr*dB_varphi_c_dr-dA_phi_dr*dB_vartheta_c_dr-d2A_phi_dr2*B_vartheta_c) &
         /bmod2-dsqg_c_dr)/sqg_c
  bder(2)=0.5d0*((dA_theta_dr*dB_varphi_c_dt-dA_phi_dr*dB_vartheta_c_dt)/bmod2-dsqg_c_dt)/sqg_c
  bder(3)=0.5d0*((dA_theta_dr*dB_varphi_c_dp-dA_phi_dr*dB_vartheta_c_dp)/bmod2-dsqg_c_dp)/sqg_c
!
  hcovar(1)=0.d0
  hcovar(2)=B_vartheta_c/bmod
  hcovar(3)=B_varphi_c/bmod
!
  hctrvr(1)=0.d0
  hctrvr(2)=Bctr_vartheta/bmod
  hctrvr(3)=Bctr_varphi/bmod
!
  hcurl(1)=((dB_varphi_c_dt-dB_vartheta_c_dp)/bmod-bder(2)*hcovar(3)+bder(3)*hcovar(2))/sqg_c
  hcurl(2)=(-dB_varphi_c_dr/bmod+bder(1)*hcovar(3))/sqg_c
  hcurl(3)=(dB_vartheta_c_dr/bmod-bder(1)*hcovar(2))/sqg_c
!
  end subroutine magfie_can
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine magfie_boozer(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
! Computes magnetic field module in units of the magnetic code  - bmod,
! square root of determinant of the metric tensor               - sqrtg,
! derivatives of the logarythm of the magnetic field module
! over coordinates                                              - bder,
! covariant componets of the unit vector of the magnetic
! field direction                                               - hcovar,
! contravariant components of this vector                       - hctrvr,
! contravariant component of the curl of this vector            - hcurl
! Order of coordinates is the following: x(1)=s (normalized toroidal flux),
! x(2)=vartheta_B (Boozer's poloidal angle), x(3)=varphi_B (Boozer's toroidal angle).
!
!  Input parameters:
!            formal:  x(3)             - array of Boozer coordinates
!  Output parameters:
!            formal:  bmod
!                     sqrtg
!                     bder(3)          - derivatives of $\log(B)$
!                     hcovar(3)        - covariant components of unit vector $\bh$ along $\bB$
!                     hctrvr(3)        - contra-variant components of unit vector $\bh$ along $\bB$
!                     hcurl(3)         - contra-variant components of curl of $\bh$
!
!  Called routines: splint_boozer_coord
!
  use vector_potentail_mod, only : torflux
!
  implicit none
!
  double precision :: bmod,sqrtg
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
!
  double precision :: r,vartheta_B,varphi_B,                                       &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,d3A_phi_dr3, &
                      B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                      B_varphi_B,dB_varphi_B,d2B_varphi_B,Bmod_B,B_r
  double precision, dimension(3) :: dBmod_B,dB_r
  double precision, dimension(6) :: d2Bmod_B,d2B_r
!
  double precision :: aiota,Bctrvr_theta,Bctrvr_phi,sqrtgbmod
!
  r=x(1)
  vartheta_B=x(2)
  varphi_B=x(3)
!
  call splint_boozer_coord(r,vartheta_B,varphi_B,                                       &
                           A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                           B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                           B_varphi_B,dB_varphi_B,d2B_varphi_B,                         &
                           Bmod_B,dBmod_B,d2Bmod_B,                                     &
                           B_r,dB_r,d2B_r)
!
  aiota=-dA_phi_dr/dA_theta_dr
!
  bmod=Bmod_B
  bder=dBmod_B/Bmod_B
!
  sqrtg=(aiota*B_vartheta_B+B_varphi_B)/bmod**2*torflux
!
  Bctrvr_phi=dA_theta_dr/sqrtg
  Bctrvr_theta=aiota*Bctrvr_phi
  hctrvr(1)=0.d0
  hctrvr(2)=Bctrvr_theta/bmod
  hctrvr(3)=Bctrvr_phi/bmod
!
  hcovar(1)=B_r/bmod
  hcovar(2)=B_vartheta_B/bmod
  hcovar(3)=B_varphi_B/bmod
!
  sqrtgbmod=sqrtg*bmod
  hcurl(1)=(B_vartheta_B*bder(3)-B_varphi_B*bder(2))/sqrtgbmod
  hcurl(2)=(B_varphi_B*bder(1)-B_r*bder(3)+dB_r(3)-dB_varphi_B)/sqrtgbmod
  hcurl(3)=(B_r*bder(2)-B_vartheta_B*bder(1)+dB_vartheta_B-dB_r(2))/sqrtgbmod
!
  end subroutine magfie_boozer
