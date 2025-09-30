module magfie_can_boozer_sub

  use libneo_kinds, only : dp
  use get_can_sub, only : splint_can_coord
  use boozer_sub, only : splint_boozer_coord
  use vector_potentail_mod, only : torflux

  implicit none

contains

  subroutine magfie_can(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
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

    implicit none

    logical :: fullset
    integer :: mode_secders
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod,sqrtg
    real(dp), intent(out) :: bder(3),hcovar(3),hctrvr(3),hcurl(3)
    real(dp) :: r,vartheta_c,varphi_c,                                       &
                        A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,d3A_phi_dr3,     &
                        sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                        B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                        B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c,     &
                        d2sqg_rr,d2sqg_rt,d2sqg_rp,d2sqg_tt,d2sqg_tp,d2sqg_pp,           &
                        d2bth_rr,d2bth_rt,d2bth_rp,d2bth_tt,d2bth_tp,d2bth_pp,           &
                        d2bph_rr,d2bph_rt,d2bph_rp,d2bph_tt,d2bph_tp,d2bph_pp
    real(dp) :: Bctr_vartheta,Bctr_varphi,bmod2

    r=x(1)
    vartheta_c=x(2)
    varphi_c=x(3)

    fullset=.false.
    mode_secders=0

    call splint_can_coord(fullset,mode_secders,r,vartheta_c,varphi_c,                      &
                          A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3,     &
                          sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                          B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                          B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                          d2sqg_rr,d2sqg_rt,d2sqg_rp,d2sqg_tt,d2sqg_tp,d2sqg_pp,           &
                          d2bth_rr,d2bth_rt,d2bth_rp,d2bth_tt,d2bth_tp,d2bth_pp,           &
                          d2bph_rr,d2bph_rt,d2bph_rp,d2bph_tt,d2bph_tp,d2bph_pp,G_c)

    sqrtg=sqg_c

    Bctr_vartheta=-dA_phi_dr/sqg_c
    Bctr_varphi=dA_theta_dr/sqg_c

    bmod2=Bctr_vartheta*B_vartheta_c+Bctr_varphi*B_varphi_c
    bmod=sqrt(bmod2)

    bder(1)=0.5d0*((dA_theta_dr*dB_varphi_c_dr-dA_phi_dr*dB_vartheta_c_dr-d2A_phi_dr2*B_vartheta_c) &
           /bmod2-dsqg_c_dr)/sqg_c
    bder(2)=0.5d0*((dA_theta_dr*dB_varphi_c_dt-dA_phi_dr*dB_vartheta_c_dt)/bmod2-dsqg_c_dt)/sqg_c
    bder(3)=0.5d0*((dA_theta_dr*dB_varphi_c_dp-dA_phi_dr*dB_vartheta_c_dp)/bmod2-dsqg_c_dp)/sqg_c

    hcovar(1)=0.d0
    hcovar(2)=B_vartheta_c/bmod
    hcovar(3)=B_varphi_c/bmod

    hctrvr(1)=0.d0
    hctrvr(2)=Bctr_vartheta/bmod
    hctrvr(3)=Bctr_varphi/bmod

    hcurl(1)=((dB_varphi_c_dt-dB_vartheta_c_dp)/bmod-bder(2)*hcovar(3)+bder(3)*hcovar(2))/sqg_c
    hcurl(2)=(-dB_varphi_c_dr/bmod+bder(1)*hcovar(3))/sqg_c
    hcurl(3)=(dB_vartheta_c_dr/bmod-bder(1)*hcovar(2))/sqg_c

  end subroutine magfie_can


  subroutine magfie_boozer(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
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

    implicit none

    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod,sqrtg
    real(dp), intent(out) :: bder(3),hcovar(3),hctrvr(3),hcurl(3)

    real(dp) :: r,vartheta_B,varphi_B,                                   &
                        A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,d3A_phi_dr3, &
                        B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                        B_varphi_B,dB_varphi_B,d2B_varphi_B,Bmod_B,B_r
    real(dp), dimension(3) :: dBmod_B,dB_r
    real(dp), dimension(6) :: d2Bmod_B,d2B_r

    real(dp) :: aiota,Bctrvr_theta,Bctrvr_phi,sqrtgbmod

    r=x(1)
    vartheta_B=x(2)
    varphi_B=x(3)

    call splint_boozer_coord(r,vartheta_B,varphi_B,                                       &
                             A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,d3A_phi_dr3, &
                             B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                             B_varphi_B,dB_varphi_B,d2B_varphi_B,                         &
                             Bmod_B,dBmod_B,d2Bmod_B,                                     &
                             B_r,dB_r,d2B_r)

    aiota=-dA_phi_dr/dA_theta_dr

    bmod=Bmod_B
    bder=dBmod_B/Bmod_B

    sqrtg=(aiota*B_vartheta_B+B_varphi_B)/bmod**2*torflux

    Bctrvr_phi=dA_theta_dr/sqrtg
    Bctrvr_theta=aiota*Bctrvr_phi
    hctrvr(1)=0.d0
    hctrvr(2)=Bctrvr_theta/bmod
    hctrvr(3)=Bctrvr_phi/bmod

    hcovar(1)=B_r/bmod
    hcovar(2)=B_vartheta_B/bmod
    hcovar(3)=B_varphi_B/bmod

    sqrtgbmod=sqrtg*bmod
    hcurl(1)=(B_vartheta_B*bder(3)-B_varphi_B*bder(2))/sqrtgbmod
    hcurl(2)=(B_varphi_B*bder(1)-B_r*bder(3)+dB_r(3)-dB_varphi_B)/sqrtgbmod
    hcurl(3)=(B_r*bder(2)-B_vartheta_B*bder(1)+dB_vartheta_B-dB_r(2))/sqrtgbmod

  end subroutine magfie_boozer

end module magfie_can_boozer_sub
