module orbit_symplectic_axis_pcart
!> Opt-in flag for crossing the magnetic axis in the symplectic Euler1 step by
!> solving its Newton in a regularized radial variable (#398).
!>
!> The flux-chart Euler1 step floors negative radial iterates to s = 0.01 and
!> reflects a converged negative radius (s,th)->(|s|,th+pi). Both corrupt the
!> implicit solve for orbits that graze the axis: the field is sampled at the
!> wrong place and the energy random-walks, so deeply trapped orbits diffuse out
!> (itpplasma/SIMPLE#398). The root cause is that the Newton Jacobian consumes the
!> field's second radial derivatives, which diverge as s^(-3/2) at the axis
!> because the m=1 harmonic scales as rho = sqrt(s); this is the s = rho^2
!> coordinate, not a field error, so no s-space healing removes it.
!>
!> Fix (planned, see #398): keep the canonical Boozer residual (so the symplectic
!> map and its conserved structure are unchanged) but solve the implicit Newton in
!> pseudo-Cartesian coordinates (X,Y) = (sqrt(s) cos theta, sqrt(s) sin theta),
!> where the field and its X,Y derivatives are finite across the axis (the sqrt(s)
!> becomes linear). The same root is found, so symplecticity is preserved; only the
!> conditioning of the solve changes. The X,Y derivatives must be evaluated
!> NATIVELY (libneo flux_pseudocartesian field layer), not chain-ruled from the
!> singular s-derivatives. Make the near-axis step fully implicit so no explicit
!> d/ds angle advance remains. Pseudo-Cartesian guiding-centre near-axis treatment:
!> Pfefferle, Cooper, Graves, Misev, Comput. Phys. Commun. 207, 144 (2016),
!> arXiv:1412.5464.

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_can_mod, only: field_can_t, eval_field => evaluate, get_derivatives
use orbit_symplectic_base, only: symplectic_integrator_t
use flux_pseudocartesian, only: pseudocart_to_flux

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

subroutine pcart_velocity(f, s, theta, phi, pphi, v)
    !> Guiding-centre velocity (Xdot, Ydot, phidot, pphidot) at (s,theta,phi,pphi),
    !> read off the canonical relations the flux-chart Euler1 uses. Finite at the
    !> axis: the 1/sqrt(s) of sdot/2rho cancels (sdot ~ s).
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: s, theta, phi, pphi
    real(dp), intent(out) :: v(4)
    real(dp) :: thdot, phidot, pphidot, pthdot, sdot, rho, cth, sth

    call eval_field(f, s, theta, phi, 0)
    call get_derivatives(f, pphi)
    thdot   = f%dH(1)/f%dpth(1)
    phidot  = (f%vpar - thdot*f%hth)/f%hph
    pphidot = -(f%dH(3) - thdot*f%dpth(3))
    pthdot  = -(f%dH(2) - thdot*f%dpth(2))
    sdot    = (pthdot - f%dpth(2)*thdot - f%dpth(3)*phidot - f%dpth(4)*pphidot) &
              / f%dpth(1)
    rho = sqrt(max(s, 0.0_dp)); cth = cos(theta); sth = sin(theta)
    if (rho > 0.0_dp) then
        v(1) = 0.5_dp*cth/rho*sdot - rho*sth*thdot
        v(2) = 0.5_dp*sth/rho*sdot + rho*cth*thdot
    else
        v(1) = 0.0_dp; v(2) = 0.0_dp
    end if
    v(3) = phidot; v(4) = pphidot
end subroutine pcart_velocity

subroutine pcart_rhs(f, y, v)
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: y(4)
    real(dp), intent(out) :: v(4)
    real(dp) :: xref(3)
    call pseudocart_to_flux([y(1), y(2), y(3)], xref)
    call pcart_velocity(f, xref(1), xref(2), y(3), y(4), v)
end subroutine pcart_rhs

subroutine axis_pcart_step(si, f, ierr)
    !> One substep of duration si%dt advanced in (X,Y,phi,pphi). The field and the
    !> guiding-centre velocity are regular across the axis in (X,Y), so the step
    !> needs no negative-s floor and no (s,theta)->(|s|,theta+pi) reflection.
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: ierr
    integer, parameter :: nsub = 8
    real(dp) :: y(4), k1(4), k2(4), k3(4), k4(4), h
    integer :: isub

    ierr = 0
    y(1) = sqrt(max(si%z(1), 0.0_dp))*cos(si%z(2))
    y(2) = sqrt(max(si%z(1), 0.0_dp))*sin(si%z(2))
    y(3) = si%z(3); y(4) = si%z(4)
    h = si%dt/real(nsub, dp)
    do isub = 1, nsub
        call pcart_rhs(f, y, k1)
        call pcart_rhs(f, y + 0.5_dp*h*k1, k2)
        call pcart_rhs(f, y + 0.5_dp*h*k2, k3)
        call pcart_rhs(f, y + h*k3, k4)
        y = y + (h/6.0_dp)*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)
    end do
    si%z(1) = y(1)**2 + y(2)**2
    si%z(2) = atan2(y(2), y(1))
    si%z(3) = y(3); si%z(4) = y(4)
    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_derivatives(f, si%z(4))
end subroutine axis_pcart_step

end module orbit_symplectic_axis_pcart
