module field_analytical_gs

use, intrinsic :: iso_fortran_env, only: dp => real64
use analytical_tokamak_field, only: &
    analytical_circular_eq_t, compute_analytical_field_cylindrical
use field_base, only: MagneticField

implicit none

type, extends(MagneticField) :: AnalyticalGSField
    type(analytical_circular_eq_t) :: eq
    real(dp) :: R0       ! Major radius [m]
    real(dp) :: a        ! Minor radius [m]
    logical :: initialized = .false.
contains
    procedure :: evaluate
end type AnalyticalGSField

contains

subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    !> Evaluate analytical GS field
    !>
    !> Note: This field uses cylindrical coordinates (R, phi, Z) instead of
    !> flux coordinates. Input x is interpreted as:
    !>   x(1) = (R - R0)/a  (normalized radial distance from axis)
    !>   x(2) = theta       (poloidal angle, unused for axisymmetric field)
    !>   x(3) = phi         (toroidal angle)
    !>
    !> Output is magnetic field in cylindrical coordinates
    class(AnalyticalGSField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    real(dp) :: R, phi, Z
    real(dp) :: B_R, B_Z, B_phi, B_mod_local
    real(dp) :: theta

    if (.not. self%initialized) then
        error stop 'AnalyticalGSField: field not initialized'
    end if

    ! Convert from normalized flux-like coordinates to cylindrical
    ! x(1) is normalized radial coordinate: rho = (R - R0)/a
    ! x(2) is poloidal angle theta
    ! x(3) is toroidal angle phi

    theta = x(2)
    phi = x(3)

    ! For circular tokamak: R = R0 + rho*a*cos(theta), Z = rho*a*sin(theta)
    R = self%R0 + x(1) * self%a * cos(theta)
    Z = x(1) * self%a * sin(theta)

    ! Evaluate field in cylindrical coordinates
    call compute_analytical_field_cylindrical(self%eq, R, phi, Z, &
                                              B_R, B_Z, B_phi, B_mod_local)

    ! For now, return field components in cylindrical basis
    ! (transformation to flux coordinates would require knowing the actual flux surfaces)
    Acov = 0.0_dp  ! Vector potential not implemented

    ! Magnetic field unit vector in cylindrical coordinates
    hcov(1) = B_R
    hcov(2) = B_phi
    hcov(3) = B_Z

    Bmod = B_mod_local

    if (present(sqgBctr)) then
        error stop 'sqgBctr not implemented for AnalyticalGSField'
    end if
end subroutine evaluate


subroutine create_analytical_gs_field(R0, epsilon, kappa, delta, A_param, B0, &
                                      Nripple, a0_ripple, alpha0, delta0, z0, &
                                      gs_field)
    !> Create analytical Grad-Shafranov field
    !>
    !> Parameters:
    !>   R0       - Major radius [m]
    !>   epsilon  - Inverse aspect ratio (a/R0)
    !>   kappa    - Elongation (optional, default 1.0)
    !>   delta    - Triangularity (optional, default 0.0)
    !>   A_param  - Shafranov parameter (optional, default -0.05)
    !>   B0       - Toroidal field on axis [T] (optional, default 1.0)
    !>   Nripple  - Number of TF coils for ripple (optional, 0 = no ripple)
    !>   a0_ripple- Ripple reference radius [m] (optional)
    !>   alpha0   - Ripple radial exponent (optional, default 2.0)
    !>   delta0   - Ripple amplitude (optional, default 0.0)
    !>   z0       - Ripple vertical center [m] (optional, default 0.0)
    real(dp), intent(in) :: R0, epsilon
    real(dp), intent(in), optional :: kappa, delta, A_param, B0
    integer, intent(in), optional :: Nripple
    real(dp), intent(in), optional :: a0_ripple, alpha0, delta0, z0
    class(AnalyticalGSField), allocatable, intent(out) :: gs_field

    real(dp) :: kappa_, delta_, A_param_, B0_
    integer :: Nripple_
    real(dp) :: a0_, alpha0_, delta0_, z0_

    ! Default values
    kappa_ = 1.0_dp
    delta_ = 0.0_dp
    A_param_ = -0.05_dp
    B0_ = 1.0_dp
    Nripple_ = 0
    a0_ = R0 * epsilon
    alpha0_ = 2.0_dp
    delta0_ = 0.0_dp
    z0_ = 0.0_dp

    ! Override with provided values
    if (present(kappa)) kappa_ = kappa
    if (present(delta)) delta_ = delta
    if (present(A_param)) A_param_ = A_param
    if (present(B0)) B0_ = B0
    if (present(Nripple)) Nripple_ = Nripple
    if (present(a0_ripple)) a0_ = a0_ripple
    if (present(alpha0)) alpha0_ = alpha0
    if (present(delta0)) delta0_ = delta0
    if (present(z0)) z0_ = z0

    allocate(AnalyticalGSField :: gs_field)

    gs_field%R0 = R0
    gs_field%a = epsilon * R0

    ! Initialize equilibrium
    if (Nripple_ > 0) then
        call gs_field%eq%init(R0, epsilon, kappa_in=kappa_, delta_in=delta_, &
                             A_in=A_param_, B0_in=B0_, &
                             Nripple_in=Nripple_, a0_in=a0_, &
                             alpha0_in=alpha0_, delta0_in=delta0_, z0_in=z0_)
    else
        call gs_field%eq%init(R0, epsilon, kappa_in=kappa_, delta_in=delta_, &
                             A_in=A_param_, B0_in=B0_)
    end if

    gs_field%initialized = .true.
end subroutine create_analytical_gs_field

end module field_analytical_gs
