module orbit_symplectic_euler1
use field_can_mod, only: field_can_t
use orbit_symplectic_base, only: symplectic_integrator_t
use, intrinsic :: iso_fortran_env, only: dp => real64

implicit none

private

public :: sympl_euler1_residual, sympl_euler1_jacobian
public :: sympl_euler1_newton_iter, sympl_euler1_extrapolate_field
public :: sympl_euler1_advance_angles

contains

subroutine sympl_euler1_residual(si, f, x, fvec)
    !$acc routine seq
    type(symplectic_integrator_t), intent(in) :: si
    type(field_can_t), intent(in) :: f
    real(dp), intent(in) :: x(2)
    real(dp), intent(out) :: fvec(2)

    fvec(1) = f%dpth(1)*(f%pth - si%pthold) &
              + si%dt*(f%dH(2)*f%dpth(1) - f%dH(1)*f%dpth(2))
    fvec(2) = f%dpth(1)*(x(2) - si%z(4)) &
              + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))
end subroutine sympl_euler1_residual

subroutine sympl_euler1_jacobian(si, f, x, jac)
    !$acc routine seq
    type(symplectic_integrator_t), intent(in) :: si
    type(field_can_t), intent(in) :: f
    real(dp), intent(in) :: x(2)
    real(dp), intent(out) :: jac(2, 2)

    jac(1, 1) = f%d2pth(1)*(f%pth - si%pthold) + f%dpth(1)**2 &
        + si%dt*(f%d2H(2)*f%dpth(1) + f%dH(2)*f%d2pth(1) &
                 - f%d2H(1)*f%dpth(2) - f%dH(1)*f%d2pth(2))
    jac(1, 2) = f%d2pth(7)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(4) &
        + si%dt*(f%d2H(8)*f%dpth(1) + f%dH(2)*f%d2pth(7) &
                 - f%d2H(7)*f%dpth(2) - f%dH(1)*f%d2pth(8))
    jac(2, 1) = f%d2pth(1)*(x(2) - si%z(4)) &
        + si%dt*(f%d2H(3)*f%dpth(1) + f%dH(3)*f%d2pth(1) &
                 - f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
    jac(2, 2) = f%d2pth(7)*(x(2) - si%z(4)) + f%dpth(1) &
        + si%dt*(f%d2H(9)*f%dpth(1) + f%dH(3)*f%d2pth(7) &
                 - f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9))
end subroutine sympl_euler1_jacobian

subroutine sympl_euler1_newton_iter(si, f, x, tolref, xlast, converged, &
        linear_failed, boundary_limited)
    !$acc routine seq
    type(symplectic_integrator_t), intent(in) :: si
    type(field_can_t), intent(in) :: f
    real(dp), intent(inout) :: x(2)
    real(dp), intent(inout) :: tolref(2)
    real(dp), intent(out) :: xlast(2)
    logical, intent(out) :: converged
    logical, intent(out) :: linear_failed
    logical, intent(out) :: boundary_limited

    real(dp) :: fvec(2), fjac(2, 2), ijac(2, 2)
    real(dp) :: correction(2), determinant, matrix_scale, step_scale

    call sympl_euler1_residual(si, f, x, fvec)
    call sympl_euler1_jacobian(si, f, x, fjac)

    determinant = fjac(1, 1)*fjac(2, 2) - fjac(1, 2)*fjac(2, 1)
    matrix_scale = maxval(abs(fjac))
    linear_failed = matrix_scale == 0d0 .or. &
        abs(determinant) <= epsilon(1d0)*matrix_scale**2
    if (linear_failed) then
        converged = .false.
        boundary_limited = .false.
        xlast = x
        return
    end if

    if (all(fjac /= 0d0)) then
        ijac(1, 1) = 1d0/(fjac(1, 1) - &
            fjac(1, 2)*fjac(2, 1)/fjac(2, 2))
        ijac(1, 2) = -1d0/(fjac(1, 1)*fjac(2, 2)/fjac(1, 2) - fjac(2, 1))
        ijac(2, 1) = -1d0/(fjac(1, 1)*fjac(2, 2)/fjac(2, 1) - fjac(1, 2))
        ijac(2, 2) = 1d0/(fjac(2, 2) - &
            fjac(1, 2)*fjac(2, 1)/fjac(1, 1))
    else
        ijac(1, 1) = fjac(2, 2)/determinant
        ijac(1, 2) = -fjac(1, 2)/determinant
        ijac(2, 1) = -fjac(2, 1)/determinant
        ijac(2, 2) = fjac(1, 1)/determinant
    end if

    xlast = x
    correction(1) = ijac(1, 1)*fvec(1) + ijac(1, 2)*fvec(2)
    correction(2) = ijac(2, 1)*fvec(1) + ijac(2, 2)*fvec(2)
    step_scale = 1d0
    if (correction(1) < 0d0) then
        step_scale = min(1d0, &
            0.8d0*max(0d0, 1d0 - x(1))/(-correction(1)))
    end if
    boundary_limited = step_scale < 1d0
    x = x - step_scale*correction

    tolref(2) = max(dabs(x(2)), tolref(2))

    converged = .not. boundary_limited .and. &
                ((dabs(fvec(1)) < si%atol .and. dabs(fvec(2)) < si%atol) .or. &
                ( &
                 dabs(x(1) - xlast(1)) < si%rtol*tolref(1) .and. &
                 dabs(x(2) - xlast(2)) < si%rtol*tolref(2)))
end subroutine sympl_euler1_newton_iter

subroutine sympl_euler1_extrapolate_field(si, f, x, xlast)
    !$acc routine seq
    type(symplectic_integrator_t), intent(in) :: si
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: x(2)
    real(dp), intent(in) :: xlast(2)

    f%pth = f%pth + f%dpth(1)*(x(1) - xlast(1)) + f%dpth(4)*(x(2) - xlast(2))
    f%dH(1) = f%dH(1) + f%d2H(1)*(x(1) - xlast(1)) + f%d2H(7)*(x(2) - xlast(2))
    f%dpth(1) = f%dpth(1) + f%d2pth(1)*(x(1) - xlast(1)) &
                + f%d2pth(7)*(x(2) - xlast(2))
    f%vpar = f%vpar + f%dvpar(1)*(x(1) - xlast(1)) + f%dvpar(4)*(x(2) - xlast(2))
    f%hth = f%hth + f%dhth(1)*(x(1) - xlast(1))
    f%hph = f%hph + f%dhph(1)*(x(1) - xlast(1))
end subroutine sympl_euler1_extrapolate_field

subroutine sympl_euler1_advance_angles(si, f)
    !$acc routine seq
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(in) :: f

    si%z(2) = si%z(2) + si%dt*f%dH(1)/f%dpth(1)
    si%z(3) = si%z(3) + si%dt*(f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph
end subroutine sympl_euler1_advance_angles

end module orbit_symplectic_euler1
