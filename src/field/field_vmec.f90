module field_vmec

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: magnetic_field_t

implicit none

type, extends(magnetic_field_t) :: vmec_field_t
contains
    procedure :: evaluate
end type vmec_field_t

contains

subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    !> Evaluate magnetic field components from VMEC equilibrium
    !> in coordinates x = (r, theta, phi) where:
    !>   r = sqrt(s) with s the normalized toroidal flux
    !>   theta = poloidal angle (NOT theta*)
    !>   phi = toroidal angle

    use spline_vmec_sub, only: splint_vmec_data, compute_field_components, metric_tensor_vmec

    class(vmec_field_t), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    real(dp) :: Acov_vartheta, Acov_varphi
    real(dp) :: dA_theta_ds, dA_phi_ds, Bctr_vartheta, Bctr_varphi
    real(dp) :: aiota, sqg, alam, dl_ds, dl_dt, dl_dp
    real(dp) :: Bcov_s, Bcov_vartheta, Bcov_varphi
    real(dp) :: s, ds_dr
    real(dp) :: gV(3, 3), sqgV

    real(dp) :: R, Z, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp

    s = x(1)**2
    ds_dr = 2d0*x(1)

    call splint_vmec_data(s, x(2), x(3), Acov_varphi, Acov_vartheta, dA_phi_ds, &
                          dA_theta_ds, aiota, R, Z, alam, dR_ds, dR_dt, dR_dp, &
                          dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)

    call compute_field_components(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                  dA_theta_ds, dA_phi_ds, dl_ds, dl_dt, dl_dp, &
                                  sqg, Bctr_vartheta, Bctr_varphi, Bcov_s, &
                                  Bcov_vartheta, Bcov_varphi)

    call metric_tensor_vmec(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, gV, sqgV)

    ! vmec_field takes (s=r**2, theta, phi) and returns symmetry flux components for A and B.
    ! Here we convert to VMEC coordinates with transformed radius, (r, theta, phi).
    Acov(1) = Acov_vartheta*dl_ds
    Acov(1) = Acov(1)*ds_dr
    Acov(2) = Acov_vartheta*(1d0 + dl_dt)
    Acov(3) = Acov_varphi + Acov_vartheta*dl_dp

    Bmod = sqrt(Bctr_vartheta*Bcov_vartheta + Bctr_varphi*Bcov_varphi)

    hcov(1) = (Bcov_s + Bcov_vartheta*dl_ds)/Bmod
    hcov(1) = hcov(1)*ds_dr
    hcov(2) = Bcov_vartheta*(1d0 + dl_dt)/Bmod
    hcov(3) = (Bcov_varphi + Bcov_vartheta*dl_dp)/Bmod

    if (present(sqgBctr)) then
        error stop "sqgBctr not implemented in vmec_field_t"
    end if
end subroutine evaluate

end module field_vmec
