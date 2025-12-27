module field_can_base

use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

real(dp), parameter :: twopi = atan(1.d0)*8.d0

integer(8) :: n_field_evaluations = 0

!$omp threadprivate(n_field_evaluations)

type :: field_can_t
    real(dp) :: Ath, Aph
    real(dp) :: hth, hph
    real(dp) :: Bmod

    real(dp), dimension(3) :: dAth, dAph
    real(dp), dimension(3) :: dhth, dhph
    real(dp), dimension(3) :: dBmod

    ! second derivatives: drdr, drdth, drdph, dthdth, dthdph, dphdph
    real(dp), dimension(6) :: d2Ath, d2Aph
    real(dp), dimension(6) :: d2hth, d2hph
    real(dp), dimension(6) :: d2Bmod

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

subroutine identity_transform(xfrom, xto)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)

    xto = xfrom
end subroutine identity_transform

end module field_can_base
