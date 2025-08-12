module magfie_sub
use spline_vmec_sub, only: vmec_field
use field_can_meiss, only: magfie_meiss
use field_can_albert, only: magfie_albert

implicit none

! Define real(dp) kind parameter
integer, parameter :: dp = kind(1.0d0)
real(dp), parameter :: twopi = 2.d0*3.14159265358979d0

abstract interface
  subroutine magfie_base(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    import :: dp
    !            x(i)   - set of 3 curvilinear space coordinates (input)
    !            bmod   - dimensionless magnetic field module: bmod=B/B_ref
    !            sqrtg  - Jacobian of space coordinates (square root of
    !                     metric tensor
    !            bder   - derivatives of logarithm of bmod over space coords
    !                     (covariant vector)
    !            hcovar - covariant components of the unit vector along
    !                     the magnetic field
    !            hctrvr - contravariant components of the unit vector along
    !                     the magnetic field
    !            hcurl  - contravariant components of the curl of this vector
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod,sqrtg
    real(dp), intent(out) :: bder(3),hcovar(3),hctrvr(3),hcurl(3)
  end subroutine magfie_base
end interface

procedure(magfie_base), pointer :: magfie => null()

integer, parameter :: TEST=-1, CANFLUX=0, VMEC=1, BOOZER=2, MEISS=3, ALBERT=4

contains

subroutine init_magfie(id)
  integer, intent(in) :: id

  select case(id)
  case(TEST)
    print *, 'init_magfie: magfie_test not implemented'
    error stop
  case(CANFLUX)
    magfie => magfie_can
  case(VMEC)
    magfie => magfie_vmec
  case(BOOZER)
    magfie => magfie_boozer
  case(MEISS)
    magfie => magfie_meiss
  case(ALBERT)
    magfie => magfie_albert
  case default
    print *,'init_magfie: unknown id ', id
    error stop
  end select
end subroutine init_magfie

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! Helper subroutine to compute finite difference derivatives
  !
  subroutine compute_vmec_derivatives(s, theta, varphi, coord_idx, step_size, &
                                     bmod_deriv, dh_ds, dh_dt, dh_dp)
    real(dp), intent(in) :: s, theta, varphi
    integer, intent(in) :: coord_idx  ! 1=s, 2=theta, 3=varphi
    real(dp), intent(in) :: step_size
    real(dp), intent(out) :: bmod_deriv
    real(dp), intent(out), optional :: dh_ds, dh_dt, dh_dp
    
    real(dp) :: s_plus, theta_plus, varphi_plus
    real(dp) :: s_minus, theta_minus, varphi_minus
    real(dp) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota
    real(dp) :: sqg, alam, dl_ds, dl_dt, dl_dp
    real(dp) :: Bctrvr_vartheta, Bctrvr_varphi
    real(dp) :: Bcovar_r, Bcovar_vartheta, Bcovar_varphi
    real(dp) :: bmod_plus, bmod_minus
    real(dp) :: bcov_s_plus, bcov_t_plus, bcov_p_plus
    real(dp) :: bcov_s_minus, bcov_t_minus, bcov_p_minus
    
    ! Initialize coordinates
    s_plus = s
    theta_plus = theta  
    varphi_plus = varphi
    s_minus = s
    theta_minus = theta
    varphi_minus = varphi
    
    ! Adjust coordinates based on which derivative we're computing
    select case(coord_idx)
    case(1)  ! s derivative
      s_plus = s + step_size
      s_minus = s - step_size
    case(2)  ! theta derivative  
      theta_plus = theta + step_size
      theta_minus = theta - step_size
    case(3)  ! varphi derivative
      varphi_plus = varphi + step_size
      varphi_minus = varphi - step_size
    end select
    
    ! Evaluate at plus step
    call vmec_field(s_plus, theta_plus, varphi_plus, A_theta, A_phi, &
                    dA_theta_ds, dA_phi_ds, aiota, sqg, alam, &
                    dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                    Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
    
    bmod_plus = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
    
    if (present(dh_ds) .or. present(dh_dt) .or. present(dh_dp)) then
      bcov_s_plus = Bcovar_r + dl_ds*Bcovar_vartheta
      bcov_t_plus = (1.d0 + dl_dt)*Bcovar_vartheta
      bcov_p_plus = Bcovar_varphi + dl_dp*Bcovar_vartheta
    end if
    
    ! Evaluate at minus step
    call vmec_field(s_minus, theta_minus, varphi_minus, A_theta, A_phi, &
                    dA_theta_ds, dA_phi_ds, aiota, sqg, alam, &
                    dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                    Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
    
    bmod_minus = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
    
    if (present(dh_ds) .or. present(dh_dt) .or. present(dh_dp)) then
      bcov_s_minus = Bcovar_r + dl_ds*Bcovar_vartheta
      bcov_t_minus = (1.d0 + dl_dt)*Bcovar_vartheta
      bcov_p_minus = Bcovar_varphi + dl_dp*Bcovar_vartheta
    end if
    
    ! Compute central differences
    bmod_deriv = (bmod_plus - bmod_minus)/(2.d0*step_size)
    
    if (present(dh_ds)) then
      dh_ds = (bcov_s_plus/bmod_plus - bcov_s_minus/bmod_minus)/(2.d0*step_size)
    end if
    
    if (present(dh_dt)) then
      dh_dt = (bcov_t_plus/bmod_plus - bcov_t_minus/bmod_minus)/(2.d0*step_size)
    end if
    
    if (present(dh_dp)) then
      dh_dp = (bcov_p_plus/bmod_plus - bcov_p_minus/bmod_minus)/(2.d0*step_size)
    end if
    
  end subroutine compute_vmec_derivatives
  
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
  !  Called routines: vmec_field, compute_vmec_derivatives
  !
  implicit none
  !
  real(dp), parameter :: hs=1.d-3, ht=hs*twopi, hp=ht/5.d0
  !
  real(dp), intent(out) :: bmod,sqrtg
  real(dp) :: s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                      sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                      Bcovar_r,Bcovar_vartheta,Bcovar_varphi
  real(dp) :: cjac,bcov_s_vmec,bcov_t_vmec,bcov_p_vmec
  real(dp) :: dhs_dt,dhs_dp,dht_ds,dht_dp,dhp_ds,dhp_dt
  real(dp), dimension(3), intent(in) :: x
  real(dp), dimension(3), intent(out) :: bder,hcovar,hctrvr,hcurl
  !
  s = x(1)
  theta = x(2)
  varphi = x(3)
  
  ! Compute derivatives using helper function
  call compute_vmec_derivatives(s, theta, varphi, 1, hs, bder(1), dht_ds, dhp_ds)
  call compute_vmec_derivatives(s, theta, varphi, 2, ht, bder(2), dhs_dt, dhp_dt)
  call compute_vmec_derivatives(s, theta, varphi, 3, hp, bder(3), dhs_dp, dht_dp)
  
  ! Evaluate at the actual point
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

  use get_can_sub, only : splint_can_coord
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
  !
  real(dp), intent(in) :: x(3)
  real(dp), intent(out) :: bmod,sqrtg
  real(dp), intent(out) :: bder(3),hcovar(3),hctrvr(3),hcurl(3)

  logical :: fullset
  integer :: mode_secders
  real(dp) :: r,vartheta_c,varphi_c,                                           &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,d3A_phi_dr3,     &
                      sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                      B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                      B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c,     &
                      d2sqg_rr,d2sqg_rt,d2sqg_rp,d2sqg_tt,d2sqg_tp,d2sqg_pp,           &
                      d2bth_rr,d2bth_rt,d2bth_rp,d2bth_tt,d2bth_tp,d2bth_pp,           &
                      d2bph_rr,d2bph_rt,d2bph_rp,d2bph_tt,d2bph_tp,d2bph_pp
  real(dp) :: Bctr_vartheta,Bctr_varphi,bmod2
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
  use boozer_sub, only : splint_boozer_coord
  !
  real(dp), intent(in) :: x(3)
  real(dp), intent(out) :: bmod,sqrtg
  real(dp), intent(out) :: bder(3),hcovar(3),hctrvr(3),hcurl(3)
  !
  real(dp) :: r,vartheta_B,varphi_B,                                       &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,d3A_phi_dr3, &
                      B_vartheta_B,dB_vartheta_B,d2B_vartheta_B,                   &
                      B_varphi_B,dB_varphi_B,d2B_varphi_B,Bmod_B,B_r
  real(dp), dimension(3) :: dBmod_B,dB_r
  real(dp), dimension(6) :: d2Bmod_B,d2B_r
  !
  real(dp) :: aiota,Bctrvr_theta,Bctrvr_phi,sqrtgbmod
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

  end module magfie_sub