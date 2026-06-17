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
use field_can_mod, only: field_can_t, eval_field => evaluate, &
  get_derivatives, get_derivatives2
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

subroutine euler1_residual_jac_q(si, f, q, pphi, theta0, r, jac)
    !> Normalized Euler1 residual r and its ANALYTIC Jacobian jac in the signed
    !> radius variables (q, pphi), s = q^2. theta is fixed at theta0 (Euler1
    !> advances it explicitly), theta0+pi for q<0 (the axis crossing). The residual
    !> is the flux-chart Euler1 residual divided by dpth/ds, which is finite at the
    !> axis. The Jacobian uses get_derivatives2 (chain/product/quotient rule, no
    !> finite differences): dr/dq = dr/ds * 2q, and dr/ds*2q is finite because
    !> dr/ds ~ s^(-1/2) while 2q ~ s^(1/2).
    type(symplectic_integrator_t), intent(in) :: si
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: q, pphi, theta0
    real(dp), intent(out) :: r(2), jac(2, 2)
    real(dp) :: s, theta, a, thd, dthd_s, dthd_p, dr1_s, dr1_p, dr2_s, dr2_p

    s = q*q
    theta = theta0
    if (q < 0.0_dp) theta = theta0 + pi
    call eval_field(f, s, theta, si%z(3), 2)
    call get_derivatives2(f, pphi)
    a = f%dpth(1)
    thd = f%dH(1)/a
    r(1) = (f%pth - si%pthold) + si%dt*(f%dH(2) - thd*f%dpth(2))
    r(2) = (pphi - si%z(4)) + si%dt*(f%dH(3) - thd*f%dpth(3))

    ! d(thd)/ds and d(thd)/dpphi by quotient rule (thd = dH(1)/dpth(1))
    dthd_s = (f%d2H(1)*a - f%dH(1)*f%d2pth(1))/(a*a)
    dthd_p = (f%d2H(7)*a - f%dH(1)*f%d2pth(7))/(a*a)
    ! dr/ds (d2 ordering: 1=ss,2=sth,3=sph,7=s pphi,8=th pphi,9=ph pphi)
    dr1_s = f%dpth(1) + si%dt*(f%d2H(2) - dthd_s*f%dpth(2) - thd*f%d2pth(2))
    dr2_s = si%dt*(f%d2H(3) - dthd_s*f%dpth(3) - thd*f%d2pth(3))
    ! dr/dpphi
    dr1_p = f%dpth(4) + si%dt*(f%d2H(8) - dthd_p*f%dpth(2) - thd*f%d2pth(8))
    dr2_p = 1.0_dp + si%dt*(f%d2H(9) - dthd_p*f%dpth(3) - thd*f%d2pth(9))
    ! chain rule s = q^2 -> column 1 scaled by 2q
    jac(1, 1) = dr1_s*2.0_dp*q; jac(1, 2) = dr1_p
    jac(2, 1) = dr2_s*2.0_dp*q; jac(2, 2) = dr2_p
end subroutine euler1_residual_jac_q

subroutine axis_pcart_step(si, f, ierr)
    !> One Euler1 substep, solved in the signed radius q (s=q^2). Same symplectic
    !> Euler1 map as the flux chart (same residual root), so there is no chart seam;
    !> the q variable just removes the negative-s floor and the (s,theta)->(|s|,
    !> theta+pi) reflection. q<0 is the axis crossing. The Newton uses the analytic
    !> Jacobian (euler1_residual_jac_q), finite at the axis. Pfefferle arXiv:1412.5464.
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: ierr
    integer, parameter :: maxit = 32
    real(dp) :: q, pphi, theta0, r(2), jac(2, 2), det, dq, dp_
    integer :: kit

    ierr = 0
    si%pthold = f%pth
    theta0 = si%z(2)
    q = sqrt(max(si%z(1), 0.0_dp))
    pphi = si%z(4)
    do kit = 1, maxit
        call euler1_residual_jac_q(si, f, q, pphi, theta0, r, jac)
        if (q*q > 1.0_dp) then
            ierr = 1
            return
        end if
        det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
        if (abs(det) < tiny(det)) then
            ierr = 1
            return
        end if
        dq = (jac(2, 2)*r(1) - jac(1, 2)*r(2))/det
        dp_ = (-jac(2, 1)*r(1) + jac(1, 1)*r(2))/det
        q = q - dq
        pphi = pphi - dp_
        if (abs(r(1)) < si%atol .and. abs(r(2)) < si%atol) exit
    end do

    si%z(1) = q*q
    si%z(4) = pphi
    if (q < 0.0_dp) si%z(2) = theta0 + pi
    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_derivatives(f, si%z(4))
    call sympl_euler1_advance_angles(si, f)
end subroutine axis_pcart_step

end module orbit_symplectic_axis_pcart
