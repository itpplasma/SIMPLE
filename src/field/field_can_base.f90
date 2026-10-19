module field_can_base

use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

real(dp), parameter :: twopi = atan(1.d0)*8.d0

integer(8) :: n_field_evaluations = 0
    real(dp) :: radial_electric_potential_slope = 0.0_dp

! Field evaluations split by derivative depth. The total above does not
! distinguish them, but the two cost very different amounts and are consumed by
! different integrator families: explicit Runge-Kutta paths only ever need first
! derivatives (mode_secders = 0), while the implicit symplectic schemes need
! second derivatives (mode_secders = 2) for their analytic Jacobians. Any
! comparison of cost per unit accuracy between those families has to weight the
! two separately, so they are counted separately.
    integer(8) :: n_field_evaluations_d1 = 0
    integer(8) :: n_field_evaluations_d2 = 0

!$omp threadprivate(n_field_evaluations)
!$omp threadprivate(n_field_evaluations_d1, n_field_evaluations_d2)

type :: field_can_t
    real(dp) :: Ath, Aph
    real(dp) :: hth, hph
    real(dp) :: Bmod
        real(dp) :: Phie = 0.0_dp

    real(dp), dimension(3) :: dAth, dAph
    real(dp), dimension(3) :: dhth, dhph
    real(dp), dimension(3) :: dBmod
        real(dp), dimension(3) :: dPhie = 0.0_dp

    ! second derivatives: drdr, drdth, drdph, dthdth, dthdph, dphdph
    real(dp), dimension(6) :: d2Ath, d2Aph
    real(dp), dimension(6) :: d2hth, d2hph
    real(dp), dimension(6) :: d2Bmod
        real(dp), dimension(6) :: d2Phie = 0.0_dp

    real(dp) :: H, pth, vpar
    real(dp), dimension(4) :: dvpar, dH, dpth

    ! order of second derivatives:
    ! d2dr2, d2drdth, d2drph, d2dth2, d2dthdph, d2dph2,
    ! d2dpphdr, d2dpphdth, d2dpphdph, d2dpph2
    real(dp), dimension(10) :: d2vpar, d2H, d2pth

    real(dp) :: mu, ro0
end type field_can_t


abstract interface
    subroutine evaluate(f, r, th_c, ph_c, mode_secders)
        import field_can_t, dp
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: r, th_c, ph_c
        integer, intent(in) :: mode_secders
    end subroutine evaluate
end interface


abstract interface
    subroutine coordinate_transform(xfrom, xto)
        import dp
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)
    end subroutine coordinate_transform
end interface

contains

    subroutine set_radial_electric_potential_slope(slope)
        real(dp), intent(in) :: slope
        radial_electric_potential_slope = slope
    end subroutine set_radial_electric_potential_slope

    subroutine apply_radial_electric_potential(f, s)
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: s

        f%Phie = radial_electric_potential_slope*s
        f%dPhie = 0.0_dp
        f%dPhie(1) = radial_electric_potential_slope
        f%d2Phie = 0.0_dp
    end subroutine apply_radial_electric_potential

! Record one field evaluation, split by derivative depth. Every field backend
! calls this instead of bumping n_field_evaluations directly, so the split can
! never drift out of step with the total.
    subroutine count_field_evaluation(mode_secders)
        integer, intent(in) :: mode_secders

        n_field_evaluations = n_field_evaluations + 1
        if (mode_secders > 0) then
            n_field_evaluations_d2 = n_field_evaluations_d2 + 1
        else
            n_field_evaluations_d1 = n_field_evaluations_d1 + 1
        end if
    end subroutine count_field_evaluation

subroutine identity_transform(xfrom, xto)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)

    xto = xfrom
end subroutine identity_transform

end module field_can_base
