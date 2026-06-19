module orbit_full_mock_cart
  ! Cartesian mock field/metric provider, no libneo dependency on the analytic
  ! paths. Flat metric (g = I), zero Christoffel. B is either uniform, a linear
  ! gradient, or (FIELD_COILS) a Biot-Savart evaluation of a coil set.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use orbit_full_provider, only: field_metric_provider_t, FO_ERR_FIELD
  use neo_biotsavart, only: coils_t, compute_magnetic_field, compute_vector_potential
  implicit none
  private
  public :: cartesian_provider_t
  public :: FIELD_UNIFORM, FIELD_LINGRAD, FIELD_COILS

  integer, parameter :: FIELD_UNIFORM = 1
  integer, parameter :: FIELD_LINGRAD = 2
  integer, parameter :: FIELD_COILS   = 3

  type, extends(field_metric_provider_t), public :: cartesian_provider_t
    integer  :: field_kind = FIELD_UNIFORM
    real(dp) :: B0(3)      = [0.0_dp, 0.0_dp, 1.0_dp]
    real(dp) :: gradB(3,3) = 0.0_dp        ! B_i(x) = B0_i + sum_j gradB(i,j) x_j
    type(coils_t), pointer :: coils => null()
  contains
    procedure :: eval_field     => cart_eval_field
    procedure :: metric         => cart_metric
    procedure :: christoffel    => cart_christoffel
    procedure :: eval_canfield  => cart_eval_canfield
    procedure :: eval_potential => cart_eval_potential
  end type cartesian_provider_t

contains

  subroutine cart_eval_field(self, x, Bvec, Bmod, hcov, ierr)
    class(cartesian_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Bvec(3), Bmod, hcov(3)
    integer,  intent(out) :: ierr
    integer :: i

    ierr = 0
    select case (self%field_kind)
    case (FIELD_UNIFORM)
      Bvec = self%B0
    case (FIELD_LINGRAD)
      do i = 1, 3
        Bvec(i) = self%B0(i) + dot_product(self%gradB(i, :), x)
      end do
    case (FIELD_COILS)
      if (.not. associated(self%coils)) then
        ierr = FO_ERR_FIELD
        Bvec = 0.0_dp; Bmod = 0.0_dp; hcov = 0.0_dp
        return
      end if
      Bvec = compute_magnetic_field(self%coils, x)
    case default
      ierr = FO_ERR_FIELD
      Bvec = 0.0_dp; Bmod = 0.0_dp; hcov = 0.0_dp
      return
    end select

    Bmod = sqrt(Bvec(1)**2 + Bvec(2)**2 + Bvec(3)**2)
    if (Bmod <= 0.0_dp) then
      ierr = FO_ERR_FIELD
      hcov = 0.0_dp
      return
    end if
    hcov = Bvec / Bmod
  end subroutine cart_eval_field

  subroutine cart_metric(self, x, g, ginv, sqrtg)
    class(cartesian_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: g(3,3), ginv(3,3), sqrtg
    integer :: i

    g = 0.0_dp
    ginv = 0.0_dp
    do i = 1, 3
      g(i,i) = 1.0_dp
      ginv(i,i) = 1.0_dp
    end do
    sqrtg = 1.0_dp
  end subroutine cart_metric

  subroutine cart_christoffel(self, x, Gamma)
    class(cartesian_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Gamma(3,3,3)

    Gamma = 0.0_dp
  end subroutine cart_christoffel

  ! Declared seam for the CPP/Pauli path; the Boris path never calls it. The
  ! mock has no analytic canonical-field model, so it flags not-implemented.
  subroutine cart_eval_canfield(self, x, Ath, Aph, hth, hph, Bmod, &
      dAth, dAph, dhth, dhph, dBmod, &
      d2Ath, d2Aph, d2hth, d2hph, d2Bmod, ierr)
    class(cartesian_provider_t), intent(in) :: self
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
  end subroutine cart_eval_canfield

  subroutine cart_eval_potential(self, x, Phi, dPhi)
    class(cartesian_provider_t), intent(in) :: self
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: Phi, dPhi(3)

    Phi = 0.0_dp
    dPhi = 0.0_dp
  end subroutine cart_eval_potential

end module orbit_full_mock_cart
