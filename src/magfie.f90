module magfie_sub
use spline_vmec_sub, only: vmec_field
use field_can_meiss, only: magfie_meiss
use field_can_albert, only: magfie_albert
use magfie_can_boozer_sub, only: magfie_can, magfie_boozer
use field_splined, only: splined_field_t
use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                              RHO_TOR, RHO_POL, PSI_TOR_NORM, PSI_POL_NORM

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

integer, parameter :: TEST=-1, CANFLUX=0, VMEC=1, BOOZER=2, MEISS=3, ALBERT=4, &
                      REFCOORDS=5

type(splined_field_t), allocatable :: refcoords_field

contains

subroutine set_magfie_refcoords_field(field)
  type(splined_field_t), intent(in) :: field

  if (allocated(refcoords_field)) deallocate (refcoords_field)
  allocate (refcoords_field, source=field)
end subroutine set_magfie_refcoords_field

subroutine init_magfie(id)
  integer, intent(in) :: id

  select case(id)
  case(TEST)
    magfie => magfie_test
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
  case(REFCOORDS)
    magfie => magfie_refcoords
  case default
    print *,'init_magfie: unknown id ', id
    error stop
  end select
end subroutine init_magfie


subroutine magfie_test(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  !> Magnetic field for analytic circular tokamak (TEST field).
  !> Coordinates: x(1)=r (minor radius), x(2)=theta (poloidal), x(3)=phi (toroidal)
  !> Uses same geometry as field_can_test: B0=1, R0=1, a=0.5, iota=1
  !>
  !> WARNING: hcurl is set to zero (curvature drift not computed).
  !> This is acceptable for symplectic integration (integmode > 0) which uses
  !> field_can_test instead. For RK45 integration (integmode=0), curvature
  !> drift would be missing - use symplectic integration with TEST field.
  implicit none

  real(dp), intent(in) :: x(3)
  real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

  real(dp), parameter :: B0 = 1.0_dp, R0 = 1.0_dp, a = 0.5_dp, iota0 = 1.0_dp
  real(dp) :: r, theta, cth, sth, R_cyl, dBmod_dr, dBmod_dth

  r = x(1)
  theta = x(2)
  cth = cos(theta)
  sth = sin(theta)

  ! Major radius at this point
  R_cyl = R0 + r * cth

  ! Magnetic field magnitude: B = B0 * (1 - r/R0 * cos(theta))
  bmod = B0 * (1.0_dp - r / R0 * cth)

  ! Jacobian sqrt(g) = r * R for circular tokamak
  sqrtg = r * R_cyl

  ! Derivatives of log(B)
  dBmod_dr = -B0 / R0 * cth
  dBmod_dth = B0 * r / R0 * sth
  bder(1) = dBmod_dr / bmod
  bder(2) = dBmod_dth / bmod
  bder(3) = 0.0_dp

  ! Covariant components of unit vector h = B/|B|
  ! In (r, theta, phi) coordinates for circular tokamak with iota=1
  hcovar(1) = 0.0_dp
  hcovar(2) = iota0 * (1.0_dp - r**2 / a**2) * r**2 / R0 / bmod
  hcovar(3) = R_cyl / bmod

  ! Contravariant components
  hctrvr(1) = 0.0_dp
  hctrvr(2) = B0 * iota0 / (r * R_cyl * bmod)
  hctrvr(3) = B0 / (r * R_cyl * bmod)

  ! Curl of h (simplified - not fully computed for TEST field)
  hcurl(1) = 0.0_dp
  hcurl(2) = 0.0_dp
  hcurl(3) = 0.0_dp

end subroutine magfie_test


subroutine magfie_refcoords(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  !> magfie in reference coordinates using analytic spline derivatives.
  !>
  !> Input x is in the splined_field_t coordinate system (typically
  !> (rho, theta, phi) where rho = sqrt(s) for VMEC coordinate systems).
  !
  implicit none

  real(dp), intent(in) :: x(3)
  real(dp), intent(out) :: bmod, sqrtg
  real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

  real(dp) :: Acov(3), Bmod_local
  real(dp) :: dAcov(3, 3), dhcov(3, 3), dBmod(3)
  real(dp) :: g(3, 3), ginv_u(3, 3), ginv_x(3, 3)
  real(dp) :: u_ref(3), J
  real(dp) :: e_cov(3, 3), jac_signed
  real(dp) :: sqrtg_abs

  if (.not. allocated(refcoords_field)) then
    print *, 'magfie_refcoords: refcoords_field not set'
    print *, 'Call set_magfie_refcoords_field before init_magfie(REFCOORDS)'
    error stop
  end if

  call refcoords_field%evaluate_with_der(x, Acov, hcovar, Bmod_local, dAcov, &
                                         dhcov, dBmod)
  bmod = Bmod_local
  bder = dBmod/max(bmod, 1.0d-30)

  call scaled_to_ref_coords(refcoords_field%coords, x, u_ref, J)
  call refcoords_field%coords%metric_tensor(u_ref, g, ginv_u, sqrtg_abs)
  call refcoords_field%coords%covariant_basis(u_ref, e_cov)

  jac_signed = signed_jacobian(e_cov)
  sqrtg = jac_signed*J

  call inverse_metric_scaled(J, ginv_u, ginv_x)
  hctrvr = matmul(ginv_x, hcovar)

  call compute_hcurl(sqrtg, dhcov, hcurl)
end subroutine magfie_refcoords


subroutine scaled_to_ref_coords(coords, x_scaled, u_ref, J)
  class(coordinate_system_t), intent(in) :: coords
  real(dp), intent(in) :: x_scaled(3)
  real(dp), intent(out) :: u_ref(3)
  real(dp), intent(out) :: J

  select type (coords)
  type is (chartmap_coordinate_system_t)
    if (coords%rho_convention == PSI_TOR_NORM .or. &
        coords%rho_convention == PSI_POL_NORM) then
      u_ref = [x_scaled(1)**2, x_scaled(2), x_scaled(3)]
      J = 2.0d0*x_scaled(1)
    else
      u_ref = x_scaled
      J = 1.0d0
    end if
  class default
    u_ref = [x_scaled(1)**2, x_scaled(2), x_scaled(3)]
    J = 2.0d0*x_scaled(1)
  end select
end subroutine scaled_to_ref_coords


subroutine inverse_metric_scaled(J, ginv_u, ginv_x)
  real(dp), intent(in) :: J
  real(dp), intent(in) :: ginv_u(3, 3)
  real(dp), intent(out) :: ginv_x(3, 3)

  ginv_x = ginv_u
  ginv_x(1, 1) = ginv_u(1, 1)/(J*J)
  ginv_x(1, 2) = ginv_u(1, 2)/J
  ginv_x(1, 3) = ginv_u(1, 3)/J
  ginv_x(2, 1) = ginv_u(2, 1)/J
  ginv_x(3, 1) = ginv_u(3, 1)/J
end subroutine inverse_metric_scaled


subroutine compute_hcurl(sqrtg, dh, hcurl)
  real(dp), intent(in) :: sqrtg
  real(dp), intent(in) :: dh(3, 3)
  real(dp), intent(out) :: hcurl(3)

  if (abs(sqrtg) <= 0.0d0) error stop 'compute_hcurl: sqrtg must be nonzero'

  hcurl(1) = (dh(2, 3) - dh(3, 2))/sqrtg
  hcurl(2) = (dh(3, 1) - dh(1, 3))/sqrtg
  hcurl(3) = (dh(1, 2) - dh(2, 1))/sqrtg
end subroutine compute_hcurl


pure function signed_jacobian(e_cov) result(jac)
  real(dp), intent(in) :: e_cov(3, 3)
  real(dp) :: jac
  real(dp) :: c(3)

  c = cross_product(e_cov(:, 2), e_cov(:, 3))
  jac = dot_product(e_cov(:, 1), c)
end function signed_jacobian


pure function cross_product(a, b) result(c)
  real(dp), intent(in) :: a(3), b(3)
  real(dp) :: c(3)

  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
end function cross_product

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

  end module magfie_sub
