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
use util, only: pi
use field_can_mod, only: field_can_t, eval_field => evaluate, get_derivatives
use orbit_symplectic_base, only: symplectic_integrator_t
use orbit_symplectic_euler1, only: sympl_euler1_advance_angles

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

subroutine euler1_residual_q(si, f, q, pphi, theta0, r)
    !> Normalized Euler1 residual at signed radius q (s=q^2). For q<0 the point is
    !> on the far side of the axis, theta0+pi. The flux-chart residual divided by
    !> dpth/ds is finite at the axis (the s^(-1/2) prefactor drops out), and only
    !> first derivatives appear, so no divergent second derivative is needed.
    type(symplectic_integrator_t), intent(in) :: si
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: q, pphi, theta0
    real(dp), intent(out) :: r(2)
    real(dp) :: s, theta, thdot

    s = q*q
    theta = theta0
    if (q < 0.0_dp) theta = theta0 + pi
    call eval_field(f, s, theta, si%z(3), 0)
    call get_derivatives(f, pphi)
    thdot = f%dH(1)/f%dpth(1)
    r(1) = (f%pth - si%pthold) + si%dt*(f%dH(2) - thdot*f%dpth(2))
    r(2) = (pphi - si%z(4)) + si%dt*(f%dH(3) - thdot*f%dpth(3))
end subroutine euler1_residual_q

subroutine axis_pcart_step(si, f, ierr)
    !> One Euler1 substep solved in the regularized variable (q, pphi), s=q^2. Same
    !> symplectic Euler1 map as the flux chart (same residual root), so there is no
    !> chart seam; only the Newton conditioning changes near the axis. The Newton
    !> uses a finite-difference Jacobian of the normalized residual, avoiding the
    !> divergent second radial derivatives the analytic Jacobian needs. A converged
    !> q<0 is the axis crossing (theta -> theta+pi); no floor, no reflection.
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: ierr
    integer, parameter :: maxit = 32
    real(dp), parameter :: hq = 1.0e-8_dp, hp = 1.0e-8_dp
    real(dp) :: q, pphi, theta0, r(2), rq(2), rp(2), jac(2, 2), det, dq, dp_
    integer :: kit
    logical :: converged

    ierr = 0
    si%pthold = f%pth          ! pth at the step start (f holds the current point)
    theta0 = si%z(2)
    q = sqrt(max(si%z(1), 0.0_dp))
    pphi = si%z(4)
    converged = .false.
    do kit = 1, maxit
        call euler1_residual_q(si, f, q, pphi, theta0, r)
        if (q*q > 1.0_dp) then
            ierr = 1
            return
        end if
        call euler1_residual_q(si, f, q + hq, pphi, theta0, rq)
        call euler1_residual_q(si, f, q, pphi + hp, theta0, rp)
        jac(:, 1) = (rq - r)/hq
        jac(:, 2) = (rp - r)/hp
        det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
        if (det == 0.0_dp) then
            ierr = 1
            return
        end if
        dq = (jac(2, 2)*r(1) - jac(1, 2)*r(2))/det
        dp_ = (-jac(2, 1)*r(1) + jac(1, 1)*r(2))/det
        q = q - dq
        pphi = pphi - dp_
        if (abs(r(1)) < si%atol .and. abs(r(2)) < si%atol) then
            converged = .true.
            exit
        end if
    end do

    si%z(1) = q*q
    si%z(4) = pphi
    if (q < 0.0_dp) si%z(2) = theta0 + pi
    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_derivatives(f, si%z(4))
    call sympl_euler1_advance_angles(si, f)
end subroutine axis_pcart_step

end module orbit_symplectic_axis_pcart
