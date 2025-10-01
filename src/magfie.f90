module magfie_sub
use spline_vmec_sub, only: vmec_field
use field_can_meiss, only: magfie_meiss
use field_can_albert, only: magfie_albert
use magfie_can_boozer_sub, only: magfie_can, magfie_boozer
use util, only: twopi
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use field_geoflux, only: geoflux_ready
use geoflux_coordinates, only: geoflux_to_cyl
use geoflux_field, only: splint_geoflux_field

implicit none

! Define real(dp) kind parameter
integer, parameter :: dp = kind(1.0d0)

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

integer, parameter :: TEST=-1, CANFLUX=0, VMEC=1, BOOZER=2, MEISS=3, ALBERT=4, GEOFLUX=5

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
    if (geoflux_ready) then
      magfie => magfie_geoflux
    else
      magfie => magfie_vmec
    end if
  case(BOOZER)
    magfie => magfie_boozer
  case(MEISS)
    magfie => magfie_meiss
  case(ALBERT)
    magfie => magfie_albert
  case(GEOFLUX)
    magfie => magfie_geoflux
  case default
    print *,'init_magfie: unknown id ', id
    error stop
  end select
end subroutine init_magfie

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
  real(dp), parameter :: twopi=2.d0*3.14159265358979d0, hs=1.d-3, ht=hs*twopi, hp=ht/5.d0
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

  subroutine magfie_geoflux(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod, sqrtg
    real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

    real(dp) :: r, theta, phi
    real(dp) :: dr_fwd, dr_bwd, dr_den
    real(dp) :: dt_step, dp_step
    real(dp) :: bmod_plus, bmod_minus
    real(dp) :: bmod_theta_plus, bmod_theta_minus
    real(dp) :: bmod_phi_plus, bmod_phi_minus
    real(dp) :: hcov_plus(3), hcov_minus(3)
    real(dp) :: hcov_theta_plus(3), hcov_theta_minus(3)
    real(dp) :: hcov_phi_plus(3), hcov_phi_minus(3)
    real(dp) :: basis(3, 3), g(3, 3), ginv(3, 3)
    real(dp) :: detg, sqrtg_geom
    real(dp) :: dh_dr(3), dh_dt(3), dh_dp(3)
    real(dp) :: phi_plus, phi_minus

    r = max(0.0_dp, min(1.0_dp, x(1)))
    theta = x(2)
    phi = x(3)

    call geoflux_eval_point(r, theta, phi, bmod, hcovar, sqrtg, basis, g, ginv, detg, sqrtg_geom)

    if (sqrtg <= 0.0_dp) sqrtg = max(sqrtg_geom, 1.0d-12)
    sqrtg = max(sqrtg, 1.0d-12)

    if (.not. ieee_is_finite(bmod)) then
      error stop 'magfie_geoflux: non-finite Bmod'
    end if
    if (.not. all(ieee_is_finite(hcovar))) then
      error stop 'magfie_geoflux: non-finite hcovar'
    end if

    dr_fwd = min(1.0d-3, 1.0_dp - r)
    dr_bwd = min(1.0d-3, r)
    dt_step = 1.0d-3*twopi
    dp_step = dt_step/5.0d0

    call geoflux_eval_basic(r + dr_fwd, theta, phi, bmod_plus, hcov_plus)
    call geoflux_eval_basic(r - dr_bwd, theta, phi, bmod_minus, hcov_minus)

    dr_den = dr_fwd + dr_bwd
    if (dr_den > 1.0d-12) then
      bder(1) = (bmod_plus - bmod_minus)/dr_den
      dh_dr = (hcov_plus - hcov_minus)/dr_den
    else
      bder(1) = 0.0_dp
      dh_dr = 0.0_dp
    end if

    call geoflux_eval_basic(r, theta + dt_step, phi, bmod_theta_plus, hcov_theta_plus)
    call geoflux_eval_basic(r, theta - dt_step, phi, bmod_theta_minus, hcov_theta_minus)
    bder(2) = (bmod_theta_plus - bmod_theta_minus)/(2.0_dp*dt_step)
    dh_dt = (hcov_theta_plus - hcov_theta_minus)/(2.0_dp*dt_step)

    phi_plus = modulo(phi + dp_step, twopi)
    phi_minus = modulo(phi - dp_step, twopi)
    call geoflux_eval_basic(r, theta, phi_plus, bmod_phi_plus, hcov_phi_plus)
    call geoflux_eval_basic(r, theta, phi_minus, bmod_phi_minus, hcov_phi_minus)
    bder(3) = (bmod_phi_plus - bmod_phi_minus)/(2.0_dp*dp_step)
    dh_dp = (hcov_phi_plus - hcov_phi_minus)/(2.0_dp*dp_step)

    bder = bder / max(bmod, 1.0d-12)

    hctrvr = matmul(ginv, hcovar)

    if (sqrtg > 0.0_dp) then
      hcurl(1) = (dh_dp(3) - dh_dt(3))/sqrtg
      hcurl(2) = (dh_dp(1) - dh_dr(3))/sqrtg
      hcurl(3) = (dh_dr(2) - dh_dt(1))/sqrtg
    else
      hcurl = 0.0_dp
    end if

  end subroutine magfie_geoflux

  subroutine geoflux_eval_point(r, theta, phi, bmod, hcov, sqrtg, basis, g, ginv, detg, sqrtg_geom)
    real(dp), intent(in) :: r, theta, phi
    real(dp), intent(out) :: bmod, hcov(3), sqrtg
    real(dp), intent(out) :: basis(3, 3), g(3, 3), ginv(3, 3)
    real(dp), intent(out) :: detg, sqrtg_geom
    real(dp) :: xcyl(3), jac(3, 3)
    real(dp) :: dRdr, dZdr, dRdtheta, dZdtheta, dRdphi, dZdphi
    real(dp) :: cosphi, sinphi, ds_dr
    real(dp) :: cross12(3)

    call geoflux_eval_basic(r, theta, phi, bmod, hcov, sqrtg, xcyl, jac)

    ds_dr = max(2.0_dp*max(r, 0.0_dp), 1.0d-8)
    cosphi = cos(xcyl(2))
    sinphi = sin(xcyl(2))

    dRdr = jac(1, 1) * ds_dr
    dZdr = jac(3, 1) * ds_dr
    dRdtheta = jac(1, 2)
    dZdtheta = jac(3, 2)
    dRdphi = jac(1, 3)
    dZdphi = jac(3, 3)

    basis(:, 1) = (/ dRdr * cosphi, dRdr * sinphi, dZdr /)
    basis(:, 2) = (/ dRdtheta * cosphi, dRdtheta * sinphi, dZdtheta /)
    basis(:, 3) = (/ dRdphi * cosphi - xcyl(1) * sinphi, &
                    dRdphi * sinphi + xcyl(1) * cosphi, dZdphi /)

    call compute_metric(basis, g, ginv, detg)
    call cross_product(basis(:, 2), basis(:, 3), cross12)
    sqrtg_geom = abs(dot_product(basis(:, 1), cross12))
    sqrtg = max(sqrtg, sqrtg_geom)
  end subroutine geoflux_eval_point

  subroutine geoflux_eval_basic(r, theta, phi, bmod, hcov, sqrtg, xcyl, jac)
    real(dp), intent(in) :: r, theta, phi
    real(dp), intent(out) :: bmod, hcov(3)
    real(dp), intent(out), optional :: sqrtg
    real(dp), intent(out), optional :: xcyl(3), jac(3, 3)

    real(dp) :: r_clip, s_geo, ds_dr
    real(dp) :: sqg_tmp(3)
    real(dp) :: acov_tmp(3), hcov_tmp(3)
    real(dp) :: cyl_tmp(3), jac_tmp(3, 3)

    r_clip = max(0.0_dp, min(1.0_dp, r))
    s_geo = r_clip * r_clip

    call geoflux_to_cyl((/ s_geo, theta, phi /), cyl_tmp, jac_tmp)
    call splint_geoflux_field(s_geo, theta, phi, acov_tmp, hcov_tmp, bmod, sqg_tmp)

    ds_dr = max(2.0_dp * max(r_clip, 0.0_dp), 1.0d-8)
    hcov(1) = hcov_tmp(1) * ds_dr
    hcov(2) = hcov_tmp(2)
    hcov(3) = hcov_tmp(3)

    if (present(sqrtg)) sqrtg = abs(sqg_tmp(1) * ds_dr)
    if (present(xcyl)) xcyl = cyl_tmp
    if (present(jac)) jac = jac_tmp
  end subroutine geoflux_eval_basic

  subroutine compute_metric(basis, g, ginv, detg)
    real(dp), intent(in) :: basis(3, 3)
    real(dp), intent(out) :: g(3, 3), ginv(3, 3)
    real(dp), intent(out) :: detg
    integer :: i, j

    do i = 1, 3
      do j = 1, 3
        g(i, j) = dot_product(basis(:, i), basis(:, j))
      end do
    end do

    call invert3x3(g, ginv, detg)
    if (abs(detg) < 1.0d-16) then
      detg = 1.0d0
      ginv = 0.0_dp
      do i = 1, 3
        ginv(i, i) = 1.0_dp
      end do
    end if
  end subroutine compute_metric

  subroutine cross_product(a, b, c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp), intent(out) :: c(3)
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end subroutine cross_product

  subroutine invert3x3(a, ainv, det)
    real(dp), intent(in) :: a(3, 3)
    real(dp), intent(out) :: ainv(3, 3)
    real(dp), intent(out) :: det

    det = a(1, 1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) &
        - a(1, 2) * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) &
        + a(1, 3) * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))

    if (abs(det) < 1.0d-16) then
      det = 0.0_dp
      ainv = 0.0_dp
      return
    end if

    ainv(1, 1) =  (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) / det
    ainv(1, 2) = -(a(1, 2) * a(3, 3) - a(1, 3) * a(3, 2)) / det
    ainv(1, 3) =  (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)) / det
    ainv(2, 1) = -(a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) / det
    ainv(2, 2) =  (a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1)) / det
    ainv(2, 3) = -(a(1, 1) * a(2, 3) - a(1, 3) * a(2, 1)) / det
    ainv(3, 1) =  (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1)) / det
    ainv(3, 2) = -(a(1, 1) * a(3, 2) - a(1, 2) * a(3, 1)) / det
    ainv(3, 3) =  (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) / det
  end subroutine invert3x3

  end module magfie_sub
