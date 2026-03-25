program test_vmec_gvec
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_vmec, only: vmec_field_t
    use field_gvec, only: gvec_field_t, create_gvec_field
    use vmec_field_eval, only: vmec_field_evaluate_with_field
    use new_vmec_stuff_mod, only: netcdffile, multharm
    use spline_vmec_sub, only: spline_vmec_data
    use params, only: pi

    implicit none

    class(vmec_field_t), allocatable :: vmec_field
    class(gvec_field_t), allocatable :: gvec_field
    real(dp) :: x(3)
    real(dp) :: Acov_vmec(3), Acov_gvec(3)
    real(dp) :: hcov_vmec(3), hcov_gvec(3)
    real(dp) :: Bmod_vmec, Bmod_gvec
    real(dp) :: max_rel_b
    real(dp) :: max_rel_a
    real(dp) :: max_rel_h
    real(dp) :: rel_err
    real(dp) :: worst_x_a(3)
    real(dp) :: worst_x_h(3)
    real(dp) :: worst_acov_vmec(3), worst_acov_gvec(3)
    real(dp) :: worst_hcov_vmec(3), worst_hcov_gvec(3)
    real(dp) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota
    real(dp) :: sqg, alam, dl_ds, dl_dt, dl_dp
    real(dp) :: Bctr_vartheta, Bctr_varphi
    real(dp) :: Bcov_s, Bcov_vartheta, Bcov_varphi
    ! The GVEC export is an interpolated interchange product, so this comparison
    ! guards convention and regression errors rather than exact spline equality.
    real(dp), parameter :: bmod_tol = 2.5e-2_dp
    real(dp), parameter :: acov_tol = 6.0e-2_dp
    real(dp), parameter :: hcov_tol = 3.0e-2_dp
    integer :: worst_a_component
    integer :: worst_h_component
    real(dp), parameter :: s_values(4) = [0.05_dp, 0.2_dp, 0.45_dp, 0.75_dp]
    real(dp), parameter :: phi_period = pi
    integer :: i, j, k

    netcdffile = 'wout.nc'
    multharm = 7
    call spline_vmec_data
    allocate (vmec_field_t :: vmec_field)
    call create_gvec_field('wout.gvec_export.nc', gvec_field)

    max_rel_b = 0.0_dp
    max_rel_a = 0.0_dp
    max_rel_h = 0.0_dp
    worst_a_component = 0
    worst_h_component = 0

    do i = 1, 4
        x(1) = s_values(i)
        do j = 1, 4
            x(2) = 2.0_dp * pi * real(j - 1, dp) / 4.0_dp
            do k = 1, 3
                x(3) = phi_period * real(k - 1, dp) / 3.0_dp

                call vmec_field%evaluate(x, Acov_vmec, hcov_vmec, Bmod_vmec)
                call gvec_field%evaluate(x, Acov_gvec, hcov_gvec, Bmod_gvec)

                rel_err = abs(Bmod_gvec - Bmod_vmec) / abs(Bmod_vmec)
                max_rel_b = max(max_rel_b, rel_err)

                rel_err = maxval(abs(Acov_gvec - Acov_vmec) / max(abs(Acov_vmec), 1.0_dp))
                if (rel_err > max_rel_a) then
                    max_rel_a = rel_err
                    worst_x_a = x
                    worst_acov_vmec = Acov_vmec
                    worst_acov_gvec = Acov_gvec
                    worst_a_component = maxloc(abs(Acov_gvec - Acov_vmec) / &
                                               max(abs(Acov_vmec), 1.0_dp), 1)
                end if

                rel_err = maxval(abs(hcov_gvec - hcov_vmec) / max(abs(hcov_vmec), 1.0_dp))
                if (rel_err > max_rel_h) then
                    max_rel_h = rel_err
                    worst_x_h = x
                    worst_hcov_vmec = hcov_vmec
                    worst_hcov_gvec = hcov_gvec
                    worst_h_component = maxloc(abs(hcov_gvec - hcov_vmec) / &
                                               max(abs(hcov_vmec), 1.0_dp), 1)
                end if
            end do
        end do
    end do

    print '(A,ES12.4)', 'max relative |B| error   = ', max_rel_b
    print '(A,ES12.4)', 'max relative Acov error  = ', max_rel_a
    print '(A,ES12.4)', 'max relative hcov error  = ', max_rel_h
    print '(A,3ES12.4)', 'worst Acov x           = ', worst_x_a
    print '(A,I0)', 'worst Acov component    = ', worst_a_component
    print '(A,3ES12.4)', 'worst Acov vmec        = ', worst_acov_vmec
    print '(A,3ES12.4)', 'worst Acov gvec        = ', worst_acov_gvec
    print '(A,3ES12.4)', 'worst hcov x           = ', worst_x_h
    print '(A,I0)', 'worst hcov component    = ', worst_h_component
    print '(A,3ES12.4)', 'worst hcov vmec        = ', worst_hcov_vmec
    print '(A,3ES12.4)', 'worst hcov gvec        = ', worst_hcov_gvec
    call vmec_field_evaluate_with_field(gvec_field, worst_x_a(1), worst_x_a(2), worst_x_a(3), &
                                        A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, sqg, &
                                        alam, dl_ds, dl_dt, dl_dp, Bctr_vartheta, &
                                        Bctr_varphi, Bcov_s, Bcov_vartheta, Bcov_varphi)
    print '(A,ES12.4)', 'worst-point dl_ds      = ', dl_ds
    print '(A,ES12.4)', 'worst-point A_theta    = ', A_theta
    print '(A,ES12.4)', 'worst-point sqg        = ', sqg

    if (max_rel_b > bmod_tol) error stop 'test_vmec_gvec: |B| mismatch'
    if (max_rel_a > acov_tol) error stop 'test_vmec_gvec: Acov mismatch'
    if (max_rel_h > hcov_tol) error stop 'test_vmec_gvec: hcov mismatch'
end program test_vmec_gvec
