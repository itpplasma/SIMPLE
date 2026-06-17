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

implicit none
private

public :: axis_pcart_enabled, axis_pcart_smax, set_axis_pcart

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

end module orbit_symplectic_axis_pcart
