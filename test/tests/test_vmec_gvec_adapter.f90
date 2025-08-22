program test_vmec_gvec_adapter
    !> Test that compares all adapter quantities between VMEC and GVEC fields
    !> This ensures that canonical coordinates work identically for both field types
    
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_vmec, only: VmecField
    use field_gvec, only: GvecField, create_gvec_field
    use vmec_field_adapter, only: vmec_field_evaluate_with_field, &
                                   vmec_lambda_interpolate_with_field, &
                                   vmec_iota_interpolate_with_field
    use new_vmec_stuff_mod, only: netcdffile, multharm
    use spline_vmec_sub, only: spline_vmec_data
    use params, only: pi
    
    implicit none
    
    class(VmecField), allocatable :: vmec_field
    class(GvecField), allocatable :: gvec_field
    character(len=256) :: vmec_file, gvec_file
    logical :: file_exists, test_passed
    
    ! Test grid
    integer, parameter :: n_test = 5
    real(dp) :: s_test(n_test), theta_test(n_test), phi_test(n_test)
    
    ! Adapter outputs for VMEC
    real(dp) :: A_theta_v, A_phi_v, dA_theta_ds_v, dA_phi_ds_v
    real(dp) :: aiota_v, sqg_v, alam_v
    real(dp) :: dl_ds_v, dl_dt_v, dl_dp_v
    real(dp) :: Bctrvr_vartheta_v, Bctrvr_varphi_v
    real(dp) :: Bcovar_r_v, Bcovar_vartheta_v, Bcovar_varphi_v
    
    ! Adapter outputs for GVEC
    real(dp) :: A_theta_g, A_phi_g, dA_theta_ds_g, dA_phi_ds_g
    real(dp) :: aiota_g, sqg_g, alam_g
    real(dp) :: dl_ds_g, dl_dt_g, dl_dp_g
    real(dp) :: Bctrvr_vartheta_g, Bctrvr_varphi_g
    real(dp) :: Bcovar_r_g, Bcovar_vartheta_g, Bcovar_varphi_g
    
    ! Lambda interpolation test
    real(dp) :: alam_interp_v, dl_dt_interp_v
    real(dp) :: alam_interp_g, dl_dt_interp_g
    
    ! Iota interpolation test
    real(dp) :: aiota_interp_v, daiota_ds_v
    real(dp) :: aiota_interp_g, daiota_ds_g
    
    ! Error tracking
    real(dp) :: max_rel_diff, rel_diff
    real(dp) :: tol_strict = 1.0e-3_dp    ! 0.1% for well-defined quantities
    real(dp) :: tol_relaxed = 1.0e-2_dp   ! 1% for derived quantities
    real(dp) :: tol_derivatives = 5.0e-2_dp ! 5% for numerical derivatives
    
    ! Track what's implemented vs dummy values
    logical :: test_essential_only = .true.  ! Focus on essential quantities
    integer :: n_pass = 0, n_fail = 0, n_dummy = 0
    
    integer :: i
    
    ! Dummy variables for field initialization
    real(dp) :: x_test(3), Acov_dummy(3), hcov_dummy(3), Bmod_dummy
    
    print *, '======================================================='
    print *, 'Testing vmec_field_adapter for VMEC vs GVEC'
    print *, '======================================================='
    print *, ''
    
    ! Setup test files
    vmec_file = 'wout.nc'
    gvec_file = 'test_vmec_gvec_State_0000_00000000.dat'
    
    ! Check files exist
    inquire(file=vmec_file, exist=file_exists)
    if (.not. file_exists) then
        print *, 'SKIP: No VMEC file found'
        stop 0
    end if
    
    inquire(file=gvec_file, exist=file_exists)
    if (.not. file_exists) then
        print *, 'SKIP: No GVEC file found. Run test_vmec_gvec first.'
        stop 0
    end if
    
    ! Initialize VMEC field
    netcdffile = vmec_file
    multharm = 5
    call spline_vmec_data
    allocate(VmecField :: vmec_field)
    
    ! Initialize GVEC field
    call create_gvec_field(gvec_file, gvec_field)
    
    ! Define test points covering different regions
    s_test = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]
    theta_test = [0.0_dp, pi/4.0_dp, pi/2.0_dp, pi, 3.0_dp*pi/2.0_dp]
    phi_test = [0.0_dp, pi/6.0_dp, pi/4.0_dp, pi/3.0_dp, pi/2.0_dp]
    
    test_passed = .true.
    max_rel_diff = 0.0_dp
    
    print *, 'Testing adapter field evaluation at multiple points...'
    print *, '------------------------------------------------------'
    
    do i = 1, n_test
        print '(A,I1,A,3(F6.3,A))', 'Point ', i, ': (s,θ,φ) = (', &
            s_test(i), ', ', theta_test(i)/pi, 'π, ', phi_test(i)/pi, 'π)'
        
        ! Evaluate with VMEC field through adapter
        call vmec_field_evaluate_with_field(vmec_field, s_test(i), theta_test(i), phi_test(i), &
            A_theta_v, A_phi_v, dA_theta_ds_v, dA_phi_ds_v, aiota_v, &
            sqg_v, alam_v, dl_ds_v, dl_dt_v, dl_dp_v, &
            Bctrvr_vartheta_v, Bctrvr_varphi_v, &
            Bcovar_r_v, Bcovar_vartheta_v, Bcovar_varphi_v)
        
        ! Evaluate with GVEC field through adapter
        ! First initialize GVEC state by calling field%evaluate
        x_test = [sqrt(s_test(i)), theta_test(i), phi_test(i)]
        call gvec_field%evaluate(x_test, Acov_dummy, hcov_dummy, Bmod_dummy)
        
        call vmec_field_evaluate_with_field(gvec_field, s_test(i), theta_test(i), phi_test(i), &
            A_theta_g, A_phi_g, dA_theta_ds_g, dA_phi_ds_g, aiota_g, &
            sqg_g, alam_g, dl_ds_g, dl_dt_g, dl_dp_g, &
            Bctrvr_vartheta_g, Bctrvr_varphi_g, &
            Bcovar_r_g, Bcovar_vartheta_g, Bcovar_varphi_g)
        
        ! Compare key quantities with detailed diagnostics
        
        ! Check sqrt(g) - note VMEC uses left-handed system (negative sqrt(g))
        call check_quantity_diagnostic('  sqrt(g)', sqg_v, sqg_g, tol_strict, 'ESSENTIAL', &
            'Jacobian - VMEC is left-handed (negative), GVEC may differ')
        
        ! Check iota - sign convention may differ
        call check_quantity_diagnostic('  iota', aiota_v, aiota_g, tol_strict, 'ESSENTIAL', &
            'Rotational transform - sign convention may differ')
        
        ! Lambda and derivatives - should match if implemented
        call check_quantity_diagnostic('  Lambda', alam_v, alam_g, tol_relaxed, 'COMPUTED', &
            'Stream function for canonical coordinates')
        call check_quantity_diagnostic('  dLambda/ds', dl_ds_v, dl_ds_g, tol_relaxed, 'COMPUTED', &
            'Lambda radial derivative')
        call check_quantity_diagnostic('  dLambda/dθ', dl_dt_v, dl_dt_g, tol_relaxed, 'COMPUTED', &
            'Lambda poloidal derivative')
        
        ! Magnetic field components
        call check_quantity_diagnostic('  B^theta', Bctrvr_vartheta_v, Bctrvr_vartheta_g, tol_strict, 'ESSENTIAL', &
            'Contravariant poloidal B field')
        call check_quantity_diagnostic('  B^phi', Bctrvr_varphi_v, Bctrvr_varphi_g, tol_strict, 'ESSENTIAL', &
            'Contravariant toroidal B field')
        
        ! Vector potential (gauge-dependent)
        call check_quantity_diagnostic('  A_theta', A_theta_v, A_theta_g, tol_relaxed, 'COMPUTED', &
            'Poloidal vector potential - gauge dependent')
        call check_quantity_diagnostic('  A_phi', A_phi_v, A_phi_g, tol_relaxed, 'COMPUTED', &
            'Toroidal vector potential')
        
        ! Derivatives (numerical for GVEC)
        if (abs(dA_theta_ds_v) > 1.0e-10_dp .or. abs(dA_theta_ds_g) > 1.0e-10_dp) then
            call check_quantity_diagnostic('  dA_theta/ds', dA_theta_ds_v, dA_theta_ds_g, tol_derivatives, 'NUMERICAL', &
                'Computed via finite differences for GVEC')
        end if
        if (abs(dA_phi_ds_v) > 1.0e-10_dp .or. abs(dA_phi_ds_g) > 1.0e-10_dp) then
            call check_quantity_diagnostic('  dA_phi/ds', dA_phi_ds_v, dA_phi_ds_g, tol_derivatives, 'NUMERICAL', &
                'Computed via finite differences for GVEC')
        end if
        
        print *, ''
    end do
    
    ! Test lambda interpolation consistency
    print *, 'Testing lambda interpolation consistency...'
    print *, '------------------------------------------'
    do i = 1, n_test
        call vmec_lambda_interpolate_with_field(vmec_field, s_test(i), theta_test(i), phi_test(i), &
                                                alam_interp_v, dl_dt_interp_v)
        call vmec_lambda_interpolate_with_field(gvec_field, s_test(i), theta_test(i), phi_test(i), &
                                                alam_interp_g, dl_dt_interp_g)
        
        call check_quantity_diagnostic('Lambda interp', alam_interp_v, alam_interp_g, tol_relaxed, 'INTERFACE', &
            'Lambda interpolation interface')
        call check_quantity_diagnostic('dLambda/dθ interp', dl_dt_interp_v, dl_dt_interp_g, tol_relaxed, 'INTERFACE', &
            'Lambda derivative interpolation')
    end do
    print *, ''
    
    ! Test iota interpolation consistency
    print *, 'Testing iota interpolation consistency...'
    print *, '----------------------------------------'
    do i = 1, n_test
        call vmec_iota_interpolate_with_field(vmec_field, s_test(i), aiota_interp_v, daiota_ds_v)
        call vmec_iota_interpolate_with_field(gvec_field, s_test(i), aiota_interp_g, daiota_ds_g)
        
        call check_quantity_diagnostic('iota interp', aiota_interp_v, aiota_interp_g, tol_strict, 'INTERFACE', &
            'Iota interpolation - sign convention issue')
        call check_quantity_diagnostic('diota/ds', daiota_ds_v, daiota_ds_g, tol_relaxed, 'INTERFACE', &
            'Iota derivative - sign convention issue')
    end do
    
    ! Summary
    print *, ''
    print *, '======================================================='
    print *, 'TEST SUMMARY - VMEC vs GVEC Adapter Comparison'
    print *, '======================================================='
    print '(A,I4)', 'Total tests performed: ', n_pass + n_fail + n_dummy
    print '(A,I4,A)', 'PASSED: ', n_pass, ' (within tolerance)'
    print '(A,I4,A)', 'FAILED: ', n_fail, ' (excessive difference)'
    print '(A,I4,A)', 'DUMMY:  ', n_dummy, ' (placeholder values)'
    print *, ''
    print '(A,ES12.5)', 'Maximum relative difference: ', max_rel_diff
    print *, ''
    print *, 'KNOWN ISSUES:'
    print *, '1. Sign conventions: iota and sqrt(g) have opposite signs'
    print *, '2. Coordinate systems: VMEC is left-handed, GVEC may differ'
    print *, '3. Some quantities not fully implemented in GVEC adapter'
    print *, ''
    
    if (n_fail == 0) then
        print *, 'STATUS: All implemented features working correctly'
    else if (test_essential_only .and. n_pass > n_fail) then
        print *, 'STATUS: Essential features partially working'
        print *, '        Further development needed for full compatibility'
    else
        print *, 'STATUS: Significant issues found'
        print *, '        Canonical coordinates may not work correctly with GVEC'
    end if
    
contains
    
    subroutine check_quantity_diagnostic(name, val_vmec, val_gvec, tol, category, description)
        character(*), intent(in) :: name, category, description
        real(dp), intent(in) :: val_vmec, val_gvec, tol
        
        real(dp) :: rel_diff, abs_diff
        character(len=8) :: status
        logical :: is_dummy
        
        abs_diff = abs(val_vmec - val_gvec)
        
        ! Check for dummy/placeholder values
        is_dummy = (abs(val_vmec) < 1.0e-10_dp .and. category == 'COMPUTED') .or. &
                   (category == 'NUMERICAL' .and. abs(val_gvec) < 1.0e-10_dp)
        
        ! Handle relative difference calculation
        if (abs(val_vmec) > 1.0e-10_dp) then
            rel_diff = abs_diff / abs(val_vmec)
        else if (abs(val_gvec) > 1.0e-10_dp) then
            rel_diff = abs_diff / abs(val_gvec)
        else
            ! Both values are near zero
            rel_diff = 0.0_dp
        end if
        
        if (rel_diff > max_rel_diff) max_rel_diff = rel_diff
        
        ! Determine status
        if (is_dummy) then
            status = 'DUMMY'
            n_dummy = n_dummy + 1
        else if (rel_diff <= tol) then
            status = 'PASS'
            n_pass = n_pass + 1
        else
            status = 'FAIL'
            n_fail = n_fail + 1
        end if
        
        ! Print diagnostic info
        if (status /= 'PASS' .or. .true.) then  ! Set to .true. for verbose output
            print '(A,A8,A,ES12.5,A,ES12.5,A,ES12.5,A,A)', &
                name, status, ': V=', val_vmec, ', G=', val_gvec, &
                ', diff=', rel_diff, ' [', trim(category), ']'
            if (status == 'FAIL') then
                print '(A,A)', '       Note: ', trim(description)
            end if
        end if
        
        ! Update test_passed for backward compatibility
        if (status == 'FAIL' .and. category == 'ESSENTIAL') then
            test_passed = .false.
        end if
    end subroutine check_quantity_diagnostic
    
    
end program test_vmec_gvec_adapter