module field_vmec

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: MagneticField
implicit none

type, extends(MagneticField) :: VmecField
contains
    procedure :: evaluate
end type VmecField

contains

subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    !> Evaluate magnetic field from VMEC equilibrium
    !> Input coordinates x = (r, theta, phi) where:
    !>   r = sqrt(s) with s the normalized toroidal flux
    !>   theta = poloidal angle (NOT theta*)
    !>   phi = toroidal angle
    !>
    !> Output field components are in (s, theta*, phi) coordinates where:
    !>   theta* = theta + Lambda(s,theta,phi)
    !>   Lambda is the stream function
    !>
    !> The transformation from (r,theta,phi) to (s,theta*,phi) components includes:
    !>   - Factor ds/dr = 2*r for the radial component
    !>   - Lambda derivatives for coupling between components
    
    use spline_vmec_sub

    class(VmecField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    real(dp) :: Acov_theta, Acov_phi
    real(dp) :: dA_theta_ds, dA_phi_ds, Bctr_theta, Bctr_phi
    real(dp) :: aiota, sqg, alam, dl_ds, dl_dt, dl_dp
    real(dp) :: Bcov_s, Bcov_theta, Bcov_phi
    real(dp) :: s, ds_dr

    s = x(1)**2
    ds_dr = 2d0*x(1)

    call vmec_field(s, x(2), x(3), Acov_theta, Acov_phi, &
        dA_theta_ds, dA_phi_ds, aiota, sqg, alam, dl_ds, dl_dt, dl_dp, &
        Bctr_theta, Bctr_phi, Bcov_s, Bcov_theta, Bcov_phi)

    Acov(1) = Acov_theta*dl_ds
    Acov(1) = Acov(1)*ds_dr
    Acov(2) = Acov_theta*(1d0 + dl_dt)
    Acov(3) = Acov_phi + Acov_theta*dl_dp

    Bmod = sqrt(Bctr_theta*Bcov_theta + Bctr_phi*Bcov_phi)

    hcov(1) = (Bcov_s + Bcov_theta*dl_ds)/Bmod
    hcov(1) = hcov(1)*ds_dr
    hcov(2) = Bcov_theta*(1d0 + dl_dt)/Bmod
    hcov(3) = (Bcov_phi + Bcov_theta*dl_dp)/Bmod

    if (present(sqgBctr)) then
        error stop "sqgBctr not implemented in VmecField"
    end if
end subroutine evaluate

end module field_vmec
