program test_field_abstraction
    !> Test that field abstraction works correctly with different field types
    !> Verifies that VMEC optimized path and generic path give same results
    
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: MagneticField
    use field, only: field_from_file
    use vmec_field_adapter, only: vmec_field_evaluate_with_field
    use params, only: pi
    
    implicit none
    
    class(MagneticField), allocatable :: vmec_field
    class(MagneticField), allocatable :: gvec_field
    character(len=256) :: vmec_file, gvec_file
    logical :: file_exists, test_passed
    
    ! Test points
    real(dp) :: s(3) = [0.1_dp, 0.5_dp, 0.9_dp]
    real(dp) :: theta(3) = [0.0_dp, pi/3.0_dp, 2.0_dp*pi/3.0_dp]
    real(dp) :: phi(3) = [0.0_dp, pi/6.0_dp, pi/4.0_dp]
    
    ! Field outputs for VMEC
    real(dp) :: A_theta_v, A_phi_v, dA_theta_ds_v, dA_phi_ds_v
    real(dp) :: aiota_v, sqg_v, alam_v
    real(dp) :: dl_ds_v, dl_dt_v, dl_dp_v
    real(dp) :: Bctrvr_vartheta_v, Bctrvr_varphi_v
    real(dp) :: Bcovar_r_v, Bcovar_vartheta_v, Bcovar_varphi_v
    
    ! Field outputs for GVEC
    real(dp) :: A_theta_g, A_phi_g, dA_theta_ds_g, dA_phi_ds_g
    real(dp) :: aiota_g, sqg_g, alam_g
    real(dp) :: dl_ds_g, dl_dt_g, dl_dp_g
    real(dp) :: Bctrvr_vartheta_g, Bctrvr_varphi_g
    real(dp) :: Bcovar_r_g, Bcovar_vartheta_g, Bcovar_varphi_g
    
    integer :: i
    real(dp) :: max_rel_diff
    
    print *, '======================================================='
    print *, 'Testing Field Abstraction with VMEC and GVEC'
    print *, '======================================================='
    print *, ''
    
    ! Check for required files
    vmec_file = 'wout.nc'
    inquire(file=vmec_file, exist=file_exists)
    if (.not. file_exists) then
        print *, 'SKIP: No VMEC file found'
        stop 0
    end if
    
    gvec_file = 'test_vmec_gvec_State_0000_00000000.dat'
    inquire(file=gvec_file, exist=file_exists)
    if (.not. file_exists) then
        print *, 'SKIP: No GVEC file found'
        stop 0
    end if
    
    ! Load fields
    vmec_field = field_from_file(trim(vmec_file))
    gvec_field = field_from_file(trim(gvec_file))
    
    test_passed = .true.
    max_rel_diff = 0.0_dp
    
    ! Test at multiple points
    do i = 1, 3
        print '(A,I1,A,3(F6.3,A))', 'Test point ', i, ': (s,θ,φ) = (', &
            s(i), ', ', theta(i), ', ', phi(i), ')'
        
        ! Evaluate with VMEC field
        call vmec_field_evaluate_with_field(vmec_field, s(i), theta(i), phi(i), &
            A_theta_v, A_phi_v, dA_theta_ds_v, dA_phi_ds_v, aiota_v, &
            sqg_v, alam_v, dl_ds_v, dl_dt_v, dl_dp_v, &
            Bctrvr_vartheta_v, Bctrvr_varphi_v, &
            Bcovar_r_v, Bcovar_vartheta_v, Bcovar_varphi_v)
        
        ! Evaluate with GVEC field
        call vmec_field_evaluate_with_field(gvec_field, s(i), theta(i), phi(i), &
            A_theta_g, A_phi_g, dA_theta_ds_g, dA_phi_ds_g, aiota_g, &
            sqg_g, alam_g, dl_ds_g, dl_dt_g, dl_dp_g, &
            Bctrvr_vartheta_g, Bctrvr_varphi_g, &
            Bcovar_r_g, Bcovar_vartheta_g, Bcovar_varphi_g)
        
        ! Compare key quantities
        call check_quantity('sqrt(g)', sqg_v, sqg_g, 1.0e-3_dp, test_passed, max_rel_diff)
        call check_quantity('iota', aiota_v, aiota_g, 1.0e-3_dp, test_passed, max_rel_diff)
        call check_quantity('Lambda', alam_v, alam_g, 1.0e-3_dp, test_passed, max_rel_diff)
        call check_quantity('B^theta', Bctrvr_vartheta_v, Bctrvr_vartheta_g, 1.0e-3_dp, test_passed, max_rel_diff)
        call check_quantity('B^phi', Bctrvr_varphi_v, Bctrvr_varphi_g, 1.0e-3_dp, test_passed, max_rel_diff)
        
        print *, ''
    end do
    
    print *, '======================================================='
    print '(A,ES12.5)', 'Maximum relative difference: ', max_rel_diff
    
    if (test_passed) then
        print *, 'TEST PASSED: Field abstraction working correctly'
        print *, 'VMEC and GVEC fields give consistent results through adapter'
    else
        print *, 'TEST FAILED: Excessive differences between VMEC and GVEC'
        error stop 1
    end if
    
contains
    
    subroutine check_quantity(name, val_vmec, val_gvec, tol, passed, max_diff)
        character(*), intent(in) :: name
        real(dp), intent(in) :: val_vmec, val_gvec, tol
        logical, intent(inout) :: passed
        real(dp), intent(inout) :: max_diff
        
        real(dp) :: rel_diff
        
        if (abs(val_vmec) > 1.0e-10_dp) then
            rel_diff = abs(val_vmec - val_gvec) / abs(val_vmec)
        else
            rel_diff = abs(val_vmec - val_gvec)
        end if
        
        if (rel_diff > max_diff) max_diff = rel_diff
        
        if (rel_diff > tol) then
            print '(A,A,A,2(ES12.5,A),ES12.5)', '  FAIL: ', name, &
                ' - VMEC: ', val_vmec, ', GVEC: ', val_gvec, &
                ', rel_diff: ', rel_diff
            passed = .false.
        else
            print '(A,A,A,ES12.5)', '  OK: ', name, ' - rel_diff: ', rel_diff
        end if
    end subroutine check_quantity
    
end program test_field_abstraction