program test_adapter_consistency
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use field_gvec, only : gvec_field_t, create_gvec_field
    use vmec_field_eval, only : vmec_data_interpolate_with_field
    use params, only : pi

    implicit none

    class(gvec_field_t), allocatable :: field
    real(dp) :: s, theta, varphi
    real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
    real(dp) :: R, Z, alam
    real(dp) :: dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
    real(dp) :: dl_ds, dl_dt, dl_dp

    call create_gvec_field('wout.gvec_export.nc', field)

    s = 0.42_dp
    theta = pi / 2.3_dp
    varphi = pi / 5.0_dp

    call vmec_data_interpolate_with_field(field, s, theta, varphi, &
                                          A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                          R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, &
                                          dZ_dp, dl_ds, dl_dt, dl_dp)

    if (A_theta == 0.0_dp) error stop 'test_adapter_consistency: A_theta is zero'
    if (dA_theta_ds == 0.0_dp) error stop 'test_adapter_consistency: dA_theta_ds is zero'
    if (abs(R) <= 0.0_dp) error stop 'test_adapter_consistency: R must be nonzero'
    if (dR_dt == 0.0_dp .or. dZ_dt == 0.0_dp) &
        error stop 'test_adapter_consistency: geometric derivatives must be nonzero'
    if (abs(dl_dt) >= 1.0_dp) error stop 'test_adapter_consistency: dLambda/dtheta unstable'

    print *, 'test_adapter_consistency passed'
end program test_adapter_consistency
