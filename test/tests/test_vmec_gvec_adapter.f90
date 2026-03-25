program test_vmec_gvec_adapter
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use field_vmec, only : vmec_field_t
    use field_gvec, only : gvec_field_t, create_gvec_field
    use vmec_field_eval, only : vmec_field_evaluate_with_field
    use new_vmec_stuff_mod, only : netcdffile, multharm
    use spline_vmec_sub, only : spline_vmec_data
    use params, only : pi

    implicit none

    class(vmec_field_t), allocatable :: vmec_field
    class(gvec_field_t), allocatable :: gvec_field
    real(dp) :: s, theta, varphi
    real(dp) :: vmec_vals(15), gvec_vals(15)
    real(dp) :: rel_err

    netcdffile = 'wout.nc'
    multharm = 7
    call spline_vmec_data
    allocate (vmec_field_t :: vmec_field)
    call create_gvec_field('wout.gvec_export.nc', gvec_field)

    s = 0.37_dp
    theta = pi / 4.0_dp
    varphi = pi / 9.0_dp

    call evaluate_point(vmec_field, s, theta, varphi, vmec_vals)
    call evaluate_point(gvec_field, s, theta, varphi, gvec_vals)

    rel_err = maxval(abs(gvec_vals - vmec_vals) / max(abs(vmec_vals), 1.0_dp))
    call print_diagnostics(vmec_vals, gvec_vals)
    print '(A,ES12.4)', 'max primitive relative error = ', rel_err

    if (rel_err > 5.0e-3_dp) error stop 'test_vmec_gvec_adapter: primitive mismatch'

contains

    subroutine print_diagnostics(vmec_vals, gvec_vals)
        real(dp), intent(in) :: vmec_vals(15)
        real(dp), intent(in) :: gvec_vals(15)

        character(len=16), parameter :: names(15) = [character(len=16) :: &
            'A_theta', 'A_phi', 'dA_theta_ds', 'dA_phi_ds', 'iota', 'sqg', &
            'Lambda', 'dLambda_ds', 'dLambda_dt', 'dLambda_dp', &
            'Bctr_vartheta', 'Bctr_varphi', 'Bcov_s', 'Bcov_vartheta', &
            'Bcov_varphi']
        real(dp) :: rels(15)
        integer :: i

        rels = abs(gvec_vals - vmec_vals) / max(abs(vmec_vals), 1.0_dp)
        do i = 1, size(names)
            print '(A16,2X,ES12.4,2X,ES12.4,2X,ES12.4)', trim(names(i)), &
                vmec_vals(i), gvec_vals(i), rels(i)
        end do
    end subroutine print_diagnostics

    subroutine evaluate_point(field, s, theta, varphi, vals)
        use field_base, only : magnetic_field_t
        class(magnetic_field_t), intent(in) :: field
        real(dp), intent(in) :: s, theta, varphi
        real(dp), intent(out) :: vals(15)

        call vmec_field_evaluate_with_field(field, s, theta, varphi, &
                                            vals(1), vals(2), vals(3), vals(4), vals(5), &
                                            vals(6), vals(7), vals(8), vals(9), vals(10), &
                                            vals(11), vals(12), vals(13), vals(14), &
                                            vals(15))
    end subroutine evaluate_point

end program test_vmec_gvec_adapter
