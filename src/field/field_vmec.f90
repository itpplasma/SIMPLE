module field_vmec
    !> VMEC equilibrium field evaluation.
    !> Uses libneo VMEC splines directly (no additional splining needed).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: magnetic_field_t
    use libneo_coordinates, only: coordinate_system_t, make_vmec_coordinate_system

    implicit none

    type, extends(magnetic_field_t) :: vmec_field_t
    contains
        procedure :: evaluate => vmec_evaluate
    end type vmec_field_t

contains

    subroutine vmec_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        !> Evaluate magnetic field from VMEC equilibrium.
        !> x = (s, theta, phi) where s = normalized toroidal flux.
        !> Returns covariant components in (s, theta, phi) coordinates.
        use spline_vmec_sub, only: splint_vmec_data, compute_field_components

        class(vmec_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out), optional :: sqgBctr(3)

        real(dp) :: Acov_vartheta, Acov_varphi
        real(dp) :: dA_theta_ds, dA_phi_ds, Bctr_vartheta, Bctr_varphi
        real(dp) :: aiota, sqg, alam, dl_ds, dl_dt, dl_dp
        real(dp) :: Bcov_s, Bcov_vartheta, Bcov_varphi
        real(dp) :: s
        real(dp) :: R, Z, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp

        s = x(1)

        call splint_vmec_data(s, x(2), x(3), Acov_varphi, Acov_vartheta, &
                              dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, &
                              dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                              dl_ds, dl_dt, dl_dp)

        call compute_field_components(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                      dA_theta_ds, dA_phi_ds, dl_ds, dl_dt, dl_dp, &
                                      sqg, Bctr_vartheta, Bctr_varphi, Bcov_s, &
                                      Bcov_vartheta, Bcov_varphi)

        Acov(1) = Acov_vartheta * dl_ds
        Acov(2) = Acov_vartheta * (1d0 + dl_dt)
        Acov(3) = Acov_varphi + Acov_vartheta * dl_dp

        Bmod = sqrt(Bctr_vartheta * Bcov_vartheta + Bctr_varphi * Bcov_varphi)

        hcov(1) = (Bcov_s + Bcov_vartheta * dl_ds) / Bmod
        hcov(2) = Bcov_vartheta * (1d0 + dl_dt) / Bmod
        hcov(3) = (Bcov_varphi + Bcov_vartheta * dl_dp) / Bmod

        if (present(sqgBctr)) then
            error stop "sqgBctr not implemented in vmec_field_t"
        end if
    end subroutine vmec_evaluate


    subroutine create_vmec_field(field)
        !> Create VMEC field with VMEC coordinate system.
        type(vmec_field_t), intent(out) :: field

        call make_vmec_coordinate_system(field%coords)
    end subroutine create_vmec_field

end module field_vmec
