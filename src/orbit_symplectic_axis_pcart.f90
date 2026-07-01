module orbit_symplectic_axis_pcart
!> Opt-in near-axis step for the symplectic integrators (#398).
!>
!> In flux coordinates s = rho^2 the m=1 harmonic scales as rho = sqrt(s), so the
!> second radial derivatives of the guiding-centre Hamiltonian diverge as
!> s^(-3/2) at the axis. The symplectic Newton consumes them, and the
!> (s,theta) -> (|s|,theta+pi) chart switch makes grazing orbits oscillate across
!> the axis; both inject non-conservative energy and spuriously lose orbits
!> (itpplasma/SIMPLE#398).
!>
!> Inside a thin shell s < smax this advances the substep with an explicit RK4 of
!> the guiding-centre equations of motion in pseudo-Cartesian coordinates
!> (X,Y) = (sqrt(s) cos theta, sqrt(s) sin theta). There the field and velocity
!> are smooth and single-valued across the axis: s = X^2 + Y^2 >= 0 never
!> reflects, so there is no floor, no chart switch, and no implicit solve. RK4 is
!> not symplectic, but the shell is thin and an orbit crosses it in a few 4th
!> order steps, so the energy error stays far below the bulk symplectic level.
!> This is the near-axis register of Pfefferle, Cooper, Graves, Misev,
!> Comput. Phys. Commun. 207, 144 (2016), arXiv:1412.5464, and matches SIMPLE's
!> RK path, which keeps these orbits confined.
!>
!> Wired only into the expl_impl_euler (integmode = 1) timestep, the default
!> and benchmark guiding-centre path; the other symplectic modes still use the
!> flux chart throughout.

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_mod, only: field_can_t, eval_field => evaluate, get_derivatives
use orbit_symplectic_base, only: symplectic_integrator_t

implicit none
private

public :: axis_pcart_enabled, axis_pcart_smax, set_axis_pcart, axis_pcart_step

logical, save :: axis_pcart_enabled = .false.
real(dp), save :: axis_pcart_smax = 0.01_dp
!$acc declare copyin(axis_pcart_enabled, axis_pcart_smax)

contains

subroutine set_axis_pcart(enabled, smax)
    logical, intent(in) :: enabled
    real(dp), intent(in), optional :: smax
    axis_pcart_enabled = enabled
    if (present(smax)) axis_pcart_smax = smax
end subroutine set_axis_pcart

subroutine pcart_velocity(f, y, fxy)
    !$acc routine seq
    !> Guiding-centre velocity in (X,Y,phi,pphi). y = (X,Y,phi,pphi). The flux
    !> velocity (sdot, thetadot, phidot, pphidot) is mapped to (Xdot, Ydot) by
    !> the transform Jacobian d(X,Y)/d(s,theta), regular for rho > 0.
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: y(4)
    real(dp), intent(out) :: fxy(4)
    real(dp) :: s, theta, rho, cth, sth, a, thd, phid, ppd, pthd, sdot

    s = y(1)*y(1) + y(2)*y(2)
    rho = sqrt(max(s, 1.0e-30_dp))
    theta = atan2(y(2), y(1))
    cth = y(1)/rho
    sth = y(2)/rho

    call eval_field(f, s, theta, y(3), 0)
    call get_derivatives(f, y(4))

    a = f%dpth(1)
    thd = f%dH(1)/a
    phid = (f%vpar - thd*f%hth)/f%hph
    ppd = -(f%dH(3) - thd*f%dpth(3))
    pthd = -(f%dH(2) - thd*f%dpth(2))
    sdot = (pthd - f%dpth(2)*thd - f%dpth(3)*phid - f%dpth(4)*ppd)/a

    fxy(1) = cth/(2.0_dp*rho)*sdot - rho*sth*thd
    fxy(2) = sth/(2.0_dp*rho)*sdot + rho*cth*thd
    fxy(3) = phid
    fxy(4) = ppd
end subroutine pcart_velocity

subroutine axis_pcart_step(si, f, ierr)
    !$acc routine seq
    !> One explicit RK4 substep in pseudo-Cartesian (X,Y,phi,pphi).
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: ierr

    real(dp) :: dt, rho0, y0(4), y(4), k1(4), k2(4), k3(4), k4(4), snew

    ierr = 0
    dt = si%dt
    rho0 = sqrt(max(si%z(1), 0.0_dp))
    y0 = [rho0*cos(si%z(2)), rho0*sin(si%z(2)), si%z(3), si%z(4)]

    call pcart_velocity(f, y0, k1)
    call pcart_velocity(f, y0 + 0.5_dp*dt*k1, k2)
    call pcart_velocity(f, y0 + 0.5_dp*dt*k2, k3)
    call pcart_velocity(f, y0 + dt*k3, k4)
    y = y0 + dt/6.0_dp*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)

    snew = y(1)*y(1) + y(2)*y(2)
    if (snew > 1.0_dp) then
        ierr = 1
        return
    end if

    si%z(1) = snew
    si%z(2) = atan2(y(2), y(1))
    si%z(3) = y(3)
    si%z(4) = y(4)
    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_derivatives(f, si%z(4))
end subroutine axis_pcart_step

end module orbit_symplectic_axis_pcart
