  module magfie_mod
! 05.03.2014
! 05.03.2014      double precision bscon,btcon,bfcon &
      double precision bxcon,bzcon,bfcon &
! 05.03.2014      ,gsr,gsf,gsz,gtr,gtf,gtz,gfr,gff,gfz
! 25.03.2016      ,gxcr,gxcf,gxcz,gzcr,gzcf,gzcz,gfr,gff,gfz
      ,gxcr,gxcfco,gxcz,gzcr,gzcfco,gzcz,gfr,gffco,gfz
! 05.03.2014 end
  end module magfie_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
! Computes magnetic field module in units of the magnetic code  - bmod,
! square root of determinant of the metric tensor               - sqrtg,
! derivatives of the logarythm of the magnetic field module
! over coordinates                                              - bder,
! covariant componets of the unit vector of the magnetic
! field direction                                               - hcovar,
! contravariant components of this vector                       - hctrvr,
! contravariant component of the curl of this vector            - hcurl
! Order of coordinates is the following: x(1)=R (big radius), 
! x(2)=phi (toroidal angle), x(3)=Z (altitude).
!
!  Input parameters:
!            formal:  x                -    array of coordinates
!  Output parameters:
!            formal:  bmod
!                     sqrtg
!                     bder
!                     hcovar
!                     hctrvr
!                     hcurl
!
!  Called routines:  GBhs,GBRZd 
!
! 05.03.2014
      use magfie_mod
! 05.03.2014 end
! 26.03.2016
      use gbpi_mod
! 26.03.2016 end
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      double precision hr,hf,hz
!
      double precision ri,fii,zi,br,bf,bz, &
      BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ, &
      BRK,BZK,BRRK,BRZK,BZRK,BZZK
! 05.03.2014
! 05.03.2014      double precision scon,tcon
      double precision xcon,zcon
!      ,bscon,btcon,bfcon &
!      ,gsr,gsf,gsz,gtr,gtf,gtz,gfr,gff,gfz
! 05.03.2014 end
! 27.03.2016
      integer ismax
      common/cismax/ismax
! 27.03.2016 end

!
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
!
! 25.03.2016
      double precision gxcf,gzcf,gff
! 25.03.2016 end
! 05.03.2014      rbig=max(x(1),1d-12)
!
!cccccc computation of gb in cylindrical co-ordinates cccccccc
! 05.03.2014 end      ri=rbig
! 05.03.2014      scon=x(1)
      xcon=x(1)
      fii=x(2)
! 05.03.2014      tcon=x(3)
      zcon=x(3)
! 05.03.2014 end      zi=x(3)

! 19.03.2010 CALL GBhs(ri,fii,zi,br,bf,bz,BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)

! 05.03.2014 CALL gbtj2(ri,fii,zi,br,bf,bz,BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)

! 05.03.2014      call gbstfp(scon,tcon,fii,bscon,btcon,bfcon &
      call gbstfp(xcon,zcon,fii,bxcon,bzcon,bfcon &
! 05.03.2014           ,gsr,gsf,gsz,gtr,gtf,gtz,gfr,gff,gfz &
           ,gxcr,gxcf,gxcz,gzcr,gzcf,gzcz,gfr,gff,gfz &
           ,ri,zi,br,bf,bz &
           ,brr,brf,brz,bfr,bff,bfz,bzr,bzf,bzz)
! 27.03.2016
      ierrfield=ismax
! 27.03.2016 end
! 26.03.2016
      if(ierrfield.eq.1) return
! 26.03.2016 end

      rbig=ri
! 25.03.2016
      gxcfco=gxcf*ri
      gzcfco=gzcf*ri
      gffco=gff*ri
! 25.03.2016 end

! 05.03.2014 end

! 18.03.2010      CALL GBRZd(ri,zi,BRK,BZK,BRRK,BRZK,BZRK,BZZK)
!
! H 14	re we add the vertical field
! 18.03.2010      br=br+BRK
! 18.03.2010      bz=bz+BZK
!
! 18.03.2010      BRR=BRR+BRRK
! 18.03.2010      BRZ=BRZ+BRZK
! 18.03.2010      BZR=BZR+BZRK
! 18.03.2010      BZZ=BZZ+BZZK
!ccccc end of gb computation cccccccccc
      bmod=dsqrt(br**2+bf**2+bz**2)
      sqrtg=rbig
      hr=br/bmod
      hf=bf/bmod
      hz=bz/bmod
!
      bder(1)=(brr*hr+bfr*hf+bzr*hz)/bmod
      bder(2)=(brf*hr+bff*hf+bzf*hz)/bmod
      bder(3)=(brz*hr+bfz*hf+bzz*hz)/bmod
!
      hcovar(1)=hr
      hcovar(2)=hf*rbig
      hcovar(3)=hz
!
      hctrvr(1)=hr
      hctrvr(2)=hf/rbig
      hctrvr(3)=hz
!
      hcurl(1)=((bzf-rbig*bfz)/bmod                        &
     +          hcovar(2)*bder(3)-hcovar(3)*bder(2))/sqrtg
      hcurl(2)=((brz-bzr)/bmod                             &
     +          hcovar(3)*bder(1)-hcovar(1)*bder(3))/sqrtg
      hcurl(3)=((bf+rbig*bfr-brf)/bmod                     &
     +          hcovar(1)*bder(2)-hcovar(2)*bder(1))/sqrtg
!
      return
      end
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
! Order of coordinates is the following: x(1)=R (big radius),
! x(2)=phi (toroidal angle), x(3)=Z (altitude).
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
  double precision :: bmod,sqrtg
  double precision :: s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                      sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                      Bcovar_r,Bcovar_vartheta,Bcovar_varphi
  double precision :: cjac,bcov_s_vmec,bcov_t_vmec,bcov_p_vmec
  double precision :: dhs_dt,dhs_dp,dht_ds,dht_dp,dhp_ds,dhp_dt
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
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
! Order of coordinates is the following: x(1)=R (big radius),
! x(2)=phi (toroidal angle), x(3)=Z (altitude).
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
!  Called routines: canonical_field
!
  implicit none
!
  logical :: fullset
  double precision :: bmod,sqrtg
  double precision :: r,vartheta_c,varphi_c,                                           &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,                 &
                      sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                      B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                      B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c
  double precision :: Bctr_vartheta,Bctr_varphi,bmod2
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
!
  r=x(1)
  vartheta_c=x(2)
  varphi_c=x(3)
!
  fullset=.false.
!
  call splint_can_coord(r,vartheta_c,varphi_c,                                           &
                        A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,                 &
                        sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                        B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                        B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                        fullset,G_c)
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
