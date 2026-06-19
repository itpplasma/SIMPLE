module orbit_full_tokamak
  ! Analytic axisymmetric circular-tokamak provider in cylindrical coordinates
  ! u = (R, phi, Z), no libneo dependency. It realizes the SAME circular-tokamak
  ! equilibrium as the field_can "test" canonical chart (field_can_test): major
  ! radius R0, minor radius a, on-axis modulus B0, constant rotational transform
  ! iota0. The full-orbit Boris pusher traces a real banana here that overlays
  ! the GC/CPP banana from the test chart.
  !
  ! Geometry: r = sqrt((R-R0)^2 + Z^2), poloidal angle th = atan2(Z, R-R0).
  ! Field in the orthonormal cylindrical frame (e_R, e_phi, e_Z):
  !   B_phi  = B0 R0 / R                       (toroidal, 1/R)
  !   B_pol  = B0 (r/R0) iota0                  along e_th = (-sin th, 0, cos th)
  ! so |B| ~ B0 R0/R to leading order, matching the test chart's
  ! Bmod = B0 (1 - r/R0 cos th) near the axis.
  !
  ! Metric is the cylindrical one (g = diag(1, R^2, 1)); the closed-form
  ! Christoffel symbols are identical to orbit_full_mock_cyl.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use orbit_full_provider, only: field_metric_provider_t, FO_ERR_FIELD
  implicit none
  private
  public :: tokamak_provider_t

  type, extends(field_metric_provider_t), public :: tokamak_provider_t
    real(dp) :: B0    = 1.0_dp
    real(dp) :: R0    = 1.0_dp
    real(dp) :: a     = 0.5_dp
    real(dp) :: iota0 = 1.0_dp
  contains
    procedure :: eval_field     => tok_eval_field
    procedure :: metric         => tok_metric
    procedure :: christoffel    => tok_christoffel
    procedure :: eval_canfield  => tok_eval_canfield
    procedure :: eval_potential => tok_eval_potential
  end type tokamak_provider_t

contains

  ! Bvec in the orthonormal cylindrical frame (e_R, e_phi, e_Z); hcov holds the
  ! covariant unit-field components h_i = B_i/|B| (h_phi^cov = R * h_phi^orth).
  subroutine tok_eval_field(self, x, Bvec, Bmod, hcov, ierr)
    class(tokamak_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Bvec(3), Bmod, hcov(3)
    integer,  intent(out) :: ierr
    real(dp) :: R, Z, dR, r_min, cth, sth, Btor, Bpol

    R = x(1); Z = x(3)
    if (R <= 0.0_dp) then
      ierr = FO_ERR_FIELD
      Bvec = 0.0_dp; Bmod = 0.0_dp; hcov = 0.0_dp
      return
    end if
    ierr = 0
    dR = R - self%R0
    r_min = sqrt(dR*dR + Z*Z)
    if (r_min > 0.0_dp) then
      cth = dR / r_min
      sth = Z / r_min
    else
      cth = 1.0_dp; sth = 0.0_dp
    end if

    Btor = self%B0 * self%R0 / R
    Bpol = self%B0 * (r_min / self%R0) * self%iota0

    ! e_th = (-sin th, 0, cos th) in (e_R, e_phi, e_Z)
    Bvec(1) = -Bpol * sth
    Bvec(2) =  Btor
    Bvec(3) =  Bpol * cth

    Bmod = sqrt(Bvec(1)**2 + Bvec(2)**2 + Bvec(3)**2)
    if (Bmod <= 0.0_dp) then
      ierr = FO_ERR_FIELD
      hcov = 0.0_dp
      return
    end if
    ! covariant unit-field: h_R = B_R/|B|, h_phi^cov = R*B_phi/|B|, h_Z = B_Z/|B|
    hcov(1) = Bvec(1) / Bmod
    hcov(2) = R * Bvec(2) / Bmod
    hcov(3) = Bvec(3) / Bmod
  end subroutine tok_eval_field

  subroutine tok_metric(self, x, g, ginv, sqrtg)
    class(tokamak_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: g(3,3), ginv(3,3), sqrtg
    real(dp) :: R

    R = x(1)
    g = 0.0_dp; ginv = 0.0_dp
    g(1,1) = 1.0_dp;  ginv(1,1) = 1.0_dp
    g(2,2) = R*R;     ginv(2,2) = 1.0_dp/(R*R)
    g(3,3) = 1.0_dp;  ginv(3,3) = 1.0_dp
    sqrtg = R
  end subroutine tok_metric

  subroutine tok_christoffel(self, x, Gamma)
    class(tokamak_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Gamma(3,3,3)
    real(dp) :: R

    R = x(1)
    Gamma = 0.0_dp
    Gamma(1,2,2) = -R
    Gamma(2,1,2) = 1.0_dp / R
    Gamma(2,2,1) = 1.0_dp / R
  end subroutine tok_christoffel

  ! CPP/Pauli seam: not exercised here (CPP runs on field_can_test).
  subroutine tok_eval_canfield(self, x, Ath, Aph, hth, hph, Bmod, &
      dAth, dAph, dhth, dhph, dBmod, &
      d2Ath, d2Aph, d2hth, d2hph, d2Bmod, ierr)
    class(tokamak_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Ath, Aph, hth, hph, Bmod
    real(dp), intent(out) :: dAth(3), dAph(3), dhth(3), dhph(3), dBmod(3)
    real(dp), intent(out) :: d2Ath(6), d2Aph(6), d2hth(6), d2hph(6), d2Bmod(6)
    integer,  intent(out) :: ierr

    Ath = 0.0_dp; Aph = 0.0_dp; hth = 0.0_dp; hph = 0.0_dp; Bmod = 0.0_dp
    dAth = 0.0_dp; dAph = 0.0_dp; dhth = 0.0_dp; dhph = 0.0_dp; dBmod = 0.0_dp
    d2Ath = 0.0_dp; d2Aph = 0.0_dp; d2hth = 0.0_dp; d2hph = 0.0_dp
    d2Bmod = 0.0_dp
    ierr = FO_ERR_FIELD
  end subroutine tok_eval_canfield

  subroutine tok_eval_potential(self, x, Phi, dPhi)
    class(tokamak_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Phi, dPhi(3)

    Phi = 0.0_dp
    dPhi = 0.0_dp
  end subroutine tok_eval_potential

end module orbit_full_tokamak
