module simple_magfie_vmec

use, intrinsic :: iso_fortran_env, only: dp => real64
use simple_magfie_base, only: MagneticField
implicit none

type, extends(MagneticField) :: VmecField
contains
    procedure :: evaluate
end type VmecField

contains

subroutine evaluate(self, x, Acov, hcov, sqgBctr, Bmod)
    use spline_vmec_sub
    class(VmecField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov, sqgBctr
    real(dp), intent(out) :: Bmod

    real(dp) :: Acov_theta, Acov_phi
    real(dp) :: dA_theta_ds, dA_phi_ds, Bctr_theta, Bctr_phi
    real(dp) :: aiota, sqg, alam, dl_ds, dl_dt, dl_dp
    real(dp) :: Bcov_s, Bcov_theta, Bcov_phi

    sqgBctr(1) = 0d0

    call vmec_field(x(1), x(2), x(3), Acov_theta, Acov_phi, &
        dA_theta_ds, dA_phi_ds, aiota, sqg, alam, dl_ds, dl_dt, dl_dp, &
        Bctr_theta, Bctr_phi, Bcov_s, Bcov_theta, Bcov_phi)

    Acov(1) = 0d0
    Acov(2) = Acov_theta
    Acov(3) = Acov_phi

    Bmod = sqrt(abs(Bctr_theta*Bcov_theta + Bctr_phi*Bcov_phi))

    hcov(1) = Bcov_s/Bmod
    hcov(2) = Bcov_theta/Bmod
    hcov(3) = Bcov_phi/Bmod

    sqgBctr(1) = 0d0
    sqgBctr(2) = -dA_phi_ds
    sqgBctr(3) = dA_theta_ds
end subroutine evaluate

end module simple_magfie_vmec
