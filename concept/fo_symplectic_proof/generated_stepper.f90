module generated_fo_stepper
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spectre_fo_gauss2, only: generated_eval => gauss2_eval, gauss2_step
    implicit none
    private
    public :: generated_eval, generated_step
contains
    subroutine generated_step(z, h, width, ro0, tol, zout, tangent, residual, ok, &
            e0, length, gmin, radial_scale)
        real(dp), intent(in) :: z(6), h, width, ro0, tol
        real(dp), intent(out) :: zout(6), tangent(6,6), residual
        logical, intent(out) :: ok
        real(dp), intent(in) :: e0, length, gmin, radial_scale
        real(dp) :: dt
        call gauss2_step(z, width, ro0, e0, length, gmin, h, zout, dt, ok, &
            tangent, radial_scale, tol, residual)
    end subroutine generated_step
end module generated_fo_stepper
