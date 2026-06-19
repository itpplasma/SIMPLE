module orbit_full_provider
  ! Abstract field/metric provider for the full-orbit pushers.
  !
  ! Kept in its own module so the Boris path carries no libneo symbol: the
  ! mocks extend this type, orbit_full only depends on the abstract interface.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: field_metric_provider_t
  public :: FO_OK, FO_ERR_FIELD, FO_ERR_NO_CONVERGE, FO_ERR_OUT_OF_DOMAIN

  ! Error codes shared by the provider seam and the integrator. Defined here so
  ! the mocks can flag a not-implemented eval_canfield without depending on
  ! orbit_full (which would create a use cycle); orbit_full re-exports them.
  integer, parameter :: FO_OK = 0
  integer, parameter :: FO_ERR_FIELD = 1        ! provider field eval failed
  integer, parameter :: FO_ERR_NO_CONVERGE = 2  ! CPP Newton did not converge
  integer, parameter :: FO_ERR_OUT_OF_DOMAIN = 3

  type, abstract :: field_metric_provider_t
  contains
    procedure(eval_field_if),     deferred :: eval_field
    procedure(metric_if),         deferred :: metric
    procedure(christoffel_if),    deferred :: christoffel
    procedure(eval_canfield_if),  deferred :: eval_canfield
    procedure(eval_potential_if), deferred :: eval_potential
  end type field_metric_provider_t

  abstract interface
    ! Physical B vector, |B|, and covariant unit-vector components h_i = B_i/|B|.
    subroutine eval_field_if(self, x, Bvec, Bmod, hcov, ierr)
      import :: field_metric_provider_t, dp
      class(field_metric_provider_t), intent(in) :: self
      real(dp), intent(in)  :: x(3)
      real(dp), intent(out) :: Bvec(3), Bmod, hcov(3)
      integer,  intent(out) :: ierr
    end subroutine eval_field_if

    ! Metric for the curvilinear push: g_ij, g^ij, sqrt(det g).
    subroutine metric_if(self, x, g, ginv, sqrtg)
      import :: field_metric_provider_t, dp
      class(field_metric_provider_t), intent(in) :: self
      real(dp), intent(in)  :: x(3)
      real(dp), intent(out) :: g(3,3), ginv(3,3), sqrtg
    end subroutine metric_if

    ! Christoffel symbols of the second kind, Gamma(i,m,n) = Gamma^i_{mn}.
    subroutine christoffel_if(self, x, Gamma)
      import :: field_metric_provider_t, dp
      class(field_metric_provider_t), intent(in) :: self
      real(dp), intent(in)  :: x(3)
      real(dp), intent(out) :: Gamma(3,3,3)
    end subroutine christoffel_if

    ! CPP / Pauli canonical field quantities, all in provider coordinates.
    ! Second-derivative layout d2?(6): order (11,12,13,22,23,33).
    subroutine eval_canfield_if(self, x, Ath, Aph, hth, hph, Bmod, &
        dAth, dAph, dhth, dhph, dBmod, &
        d2Ath, d2Aph, d2hth, d2hph, d2Bmod, ierr)
      import :: field_metric_provider_t, dp
      class(field_metric_provider_t), intent(in) :: self
      real(dp), intent(in)  :: x(3)
      real(dp), intent(out) :: Ath, Aph, hth, hph, Bmod
      real(dp), intent(out) :: dAth(3), dAph(3), dhth(3), dhph(3), dBmod(3)
      real(dp), intent(out) :: d2Ath(6), d2Aph(6), d2hth(6), d2hph(6), d2Bmod(6)
      integer,  intent(out) :: ierr
    end subroutine eval_canfield_if

    ! Electrostatic potential Phi and covariant gradient. Mocks return 0.
    subroutine eval_potential_if(self, x, Phi, dPhi)
      import :: field_metric_provider_t, dp
      class(field_metric_provider_t), intent(in) :: self
      real(dp), intent(in)  :: x(3)
      real(dp), intent(out) :: Phi, dPhi(3)
    end subroutine eval_potential_if
  end interface
end module orbit_full_provider
