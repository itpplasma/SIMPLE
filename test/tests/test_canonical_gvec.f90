program test_canonical_gvec
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use field_gvec, only : gvec_field_t, create_gvec_field
    use vmec_field_eval, only : vmec_field_evaluate_with_field, &
                                vmec_lambda_interpolate_with_field, &
                                vmec_iota_interpolate_with_field
    use params, only : pi

    implicit none

    class(gvec_field_t), allocatable :: field
    real(dp) :: x(3)
    real(dp) :: Acov(3), hcov(3), Bmod
    real(dp) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota
    real(dp) :: sqg, alam, dl_ds, dl_dt, dl_dp
    real(dp) :: Bctr_vartheta, Bctr_varphi
    real(dp) :: Bcov_s, Bcov_vartheta, Bcov_varphi
    real(dp) :: alam_interp, dl_dt_interp
    real(dp) :: aiota_interp, daiota_ds
    real(dp) :: Acov_from_vmec(3), hcov_from_vmec(3), Bmod_from_vmec

    call create_gvec_field('wout.gvec_export.nc', field)

    x = [0.31_dp, pi / 3.0_dp, pi / 7.0_dp]
    call field%evaluate(x, Acov, hcov, Bmod)

    call vmec_field_evaluate_with_field(field, x(1), x(2), x(3), &
                                        A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                                        sqg, alam, dl_ds, dl_dt, dl_dp, &
                                        Bctr_vartheta, Bctr_varphi, &
                                        Bcov_s, Bcov_vartheta, Bcov_varphi)

    Acov_from_vmec(1) = A_theta * dl_ds
    Acov_from_vmec(2) = A_theta * (1.0_dp + dl_dt)
    Acov_from_vmec(3) = A_phi + A_theta * dl_dp
    Bmod_from_vmec = sqrt(Bctr_vartheta * Bcov_vartheta + Bctr_varphi * Bcov_varphi)
    hcov_from_vmec(1) = (Bcov_s + Bcov_vartheta * dl_ds) / Bmod_from_vmec
    hcov_from_vmec(2) = Bcov_vartheta * (1.0_dp + dl_dt) / Bmod_from_vmec
    hcov_from_vmec(3) = (Bcov_varphi + Bcov_vartheta * dl_dp) / Bmod_from_vmec

    call vmec_lambda_interpolate_with_field(field, x(1), x(2), x(3), alam_interp, &
                                            dl_dt_interp)
    call vmec_iota_interpolate_with_field(field, x(1), aiota_interp, daiota_ds)

    if (maxval(abs(Acov - Acov_from_vmec)) > 1.0e-8_dp * max(1.0_dp, maxval(abs(Acov)))) &
        error stop 'test_canonical_gvec: Acov reconstruction mismatch'
    if (maxval(abs(hcov - hcov_from_vmec)) > 1.0e-8_dp * max(1.0_dp, maxval(abs(hcov)))) &
        error stop 'test_canonical_gvec: hcov reconstruction mismatch'
    if (abs(Bmod - Bmod_from_vmec) > 1.0e-8_dp * max(1.0_dp, abs(Bmod))) &
        error stop 'test_canonical_gvec: Bmod reconstruction mismatch'
    if (abs(alam - alam_interp) > 1.0e-10_dp) &
        error stop 'test_canonical_gvec: Lambda interpolation mismatch'
    if (abs(dl_dt - dl_dt_interp) > 1.0e-10_dp) &
        error stop 'test_canonical_gvec: dLambda/dtheta mismatch'
    if (abs(aiota - aiota_interp) > 1.0e-10_dp) &
        error stop 'test_canonical_gvec: iota interpolation mismatch'
    if (abs(daiota_ds) <= 0.0_dp) &
        error stop 'test_canonical_gvec: iota derivative should be nonzero'

    print *, 'test_canonical_gvec passed'
end program test_canonical_gvec
