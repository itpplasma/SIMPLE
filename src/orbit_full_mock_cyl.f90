module orbit_full_mock_cyl
  ! Cylindrical mock provider in coordinates u = (R, phi, Z), no libneo
  ! dependency. Analytic metric and closed-form Christoffel symbols, plus an
  ! axisymmetric purely-toroidal 1/R field for the curvature + grad-B drift
  ! oracle.
  !
  !   g    = diag(1, R^2, 1),  ginv = diag(1, 1/R^2, 1),  sqrtg = R
  !   B    = B0*R0/R in the orthonormal toroidal direction (|B| = B0*R0/R)
  !   B_phi^cov = |B|*R = B0*R0  (covariant toroidal component)
  !
  ! Nonzero Christoffel (second kind), Gamma(i,m,n) = Gamma^i_{mn}:
  !   Gamma^R_{phi phi}  = Gamma(1,2,2) = -R
  !   Gamma^phi_{R phi}  = Gamma(2,1,2) =  1/R
  !   Gamma^phi_{phi R}  = Gamma(2,2,1) =  1/R
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use orbit_full_provider, only: field_metric_provider_t, FO_ERR_FIELD
  implicit none
  private
  public :: cylindrical_provider_t

  type, extends(field_metric_provider_t), public :: cylindrical_provider_t
    real(dp) :: B0 = 1.0_dp
    real(dp) :: R0 = 1.0_dp
  contains
    procedure :: eval_field     => cyl_eval_field
    procedure :: metric         => cyl_metric
    procedure :: christoffel    => cyl_christoffel
    procedure :: eval_canfield  => cyl_eval_canfield
    procedure :: eval_potential => cyl_eval_potential
  end type cylindrical_provider_t

contains

  ! Bvec returned in the orthonormal cylindrical frame (e_R, e_phi, e_Z);
  ! hcov holds the covariant unit-field components h_i = B_i/|B|, so that
  ! h_phi = R for this purely-toroidal field.
  subroutine cyl_eval_field(self, x, Bvec, Bmod, hcov, ierr)
    class(cylindrical_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Bvec(3), Bmod, hcov(3)
    integer,  intent(out) :: ierr
    real(dp) :: R

    R = x(1)
    if (R <= 0.0_dp) then
      ierr = FO_ERR_FIELD
      Bvec = 0.0_dp; Bmod = 0.0_dp; hcov = 0.0_dp
      return
    end if
    ierr = 0
    Bmod = self%B0 * self%R0 / R
    Bvec = [0.0_dp, Bmod, 0.0_dp]          ! orthonormal toroidal
    hcov = [0.0_dp, R, 0.0_dp]             ! h_phi^cov = B_phi^cov/|B| = R
  end subroutine cyl_eval_field

  subroutine cyl_metric(self, x, g, ginv, sqrtg)
    class(cylindrical_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: g(3,3), ginv(3,3), sqrtg
    real(dp) :: R

    R = x(1)
    g = 0.0_dp
    ginv = 0.0_dp
    g(1,1) = 1.0_dp;     ginv(1,1) = 1.0_dp
    g(2,2) = R*R;        ginv(2,2) = 1.0_dp / (R*R)
    g(3,3) = 1.0_dp;     ginv(3,3) = 1.0_dp
    sqrtg = R
  end subroutine cyl_metric

  subroutine cyl_christoffel(self, x, Gamma)
    class(cylindrical_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Gamma(3,3,3)
    real(dp) :: R

    R = x(1)
    Gamma = 0.0_dp
    Gamma(1,2,2) = -R
    Gamma(2,1,2) = 1.0_dp / R
    Gamma(2,2,1) = 1.0_dp / R
  end subroutine cyl_christoffel

  ! Declared CPP/Pauli seam; not exercised by the Boris path. Flags
  ! not-implemented rather than returning a partial canonical model.
  subroutine cyl_eval_canfield(self, x, Ath, Aph, hth, hph, Bmod, &
      dAth, dAph, dhth, dhph, dBmod, &
      d2Ath, d2Aph, d2hth, d2hph, d2Bmod, ierr)
    class(cylindrical_provider_t), intent(in) :: self
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
  end subroutine cyl_eval_canfield

  subroutine cyl_eval_potential(self, x, Phi, dPhi)
    class(cylindrical_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Phi, dPhi(3)

    Phi = 0.0_dp
    dPhi = 0.0_dp
  end subroutine cyl_eval_potential

end module orbit_full_mock_cyl
