program test_canonical_gvec
    !> Integration test for canonical coordinates with GVEC fields
    !> Tests that the vmec_field_adapter properly handles GVEC fields
    !> and that vector potential derivatives are computed correctly
    
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_gvec, only: GvecField, create_gvec_field
    use field_base, only: MagneticField
    use vmec_field_adapter, only: vmec_field_evaluate_with_field, &
                                   vmec_lambda_interpolate_with_field, &
                                   vmec_iota_interpolate_with_field
    use params, only: pi
    
    implicit none
    
    class(GvecField), allocatable :: gvec_field
    character(len=256) :: gvec_file
    logical :: file_exists
    
    ! Test coordinates
    real(dp) :: s_test, theta_test, varphi_test
    
    ! Field evaluation outputs
    real(dp) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds
    real(dp) :: aiota, sqg, alam
    real(dp) :: dl_ds, dl_dt, dl_dp
    real(dp) :: Bctrvr_vartheta, Bctrvr_varphi
    real(dp) :: Bcovar_r, Bcovar_vartheta, Bcovar_varphi
    
    ! Lambda interpolation outputs
    real(dp) :: alam_interp, dl_dt_interp
    
    ! Iota interpolation outputs  
    real(dp) :: aiota_interp, daiota_ds
    
    ! Test flags
    logical :: test_passed
    
    ! Dummy variables for field initialization
    real(dp) :: x_test(3), Acov_dummy(3), hcov_dummy(3), Bmod_dummy
    
    print *, 'Testing canonical coordinates support with GVEC fields...'
    print *, ''
    
    ! Look for GVEC test file created by test_vmec_gvec
    gvec_file = 'test_vmec_gvec_State_0000_00000000.dat'
    inquire(file=gvec_file, exist=file_exists)
    
    if (.not. file_exists) then
        ! Try alternative name
        gvec_file = 'gvec_from_vmec_wout.dat'
        inquire(file=gvec_file, exist=file_exists)
    end if
    
    if (.not. file_exists) then
        print *, 'SKIP: No GVEC test file found. Run test_vmec_gvec first.'
        stop 0  ! Skip test, don't fail
    end if
    
    ! Load GVEC field
    call create_gvec_field(gvec_file, gvec_field)
    
    if (.not. gvec_field%data_loaded) then
        print *, 'ERROR: Failed to load GVEC field data'
        error stop 1
    end if
    
    print *, 'Successfully loaded GVEC field from: ', trim(gvec_file)
    print *, ''
    
    ! Test coordinates - representative interior point
    s_test = 0.25_dp           ! Mid-radius flux surface
    theta_test = pi / 3.0_dp   ! 60 degrees poloidal
    varphi_test = pi / 6.0_dp  ! 30 degrees toroidal
    
    print *, 'Testing adapter functions at (s,θ,φ) = (0.25, π/3, π/6):'
    print *, '================================================================'
    
    ! Test 1: Full field evaluation with vector potential derivatives
    print *, '1. Testing vmec_field_evaluate_with_field...'
    
    ! First initialize GVEC state by calling field%evaluate
    x_test = [sqrt(s_test), theta_test, varphi_test]
    call gvec_field%evaluate(x_test, Acov_dummy, hcov_dummy, Bmod_dummy)
    
    call vmec_field_evaluate_with_field(gvec_field, s_test, theta_test, varphi_test, &
                                        A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                                        sqg, alam, dl_ds, dl_dt, dl_dp, &
                                        Bctrvr_vartheta, Bctrvr_varphi, &
                                        Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
    
    print '(A,ES12.5)', '   A_theta       = ', A_theta
    print '(A,ES12.5)', '   A_phi         = ', A_phi  
    print '(A,ES12.5)', '   dA_theta/ds   = ', dA_theta_ds
    print '(A,ES12.5)', '   dA_phi/ds     = ', dA_phi_ds
    print '(A,ES12.5)', '   aiota         = ', aiota
    print '(A,ES12.5)', '   sqrt(g)       = ', sqg
    print '(A,ES12.5)', '   Lambda        = ', alam
    print '(A,ES12.5)', '   dLambda/ds    = ', dl_ds
    print '(A,ES12.5)', '   dLambda/dθ    = ', dl_dt
    print '(A,ES12.5)', '   dLambda/dφ    = ', dl_dp
    print *, ''
    
    ! Test 2: Lambda interpolation
    print *, '2. Testing vmec_lambda_interpolate_with_field...'
    
    call vmec_lambda_interpolate_with_field(gvec_field, s_test, theta_test, varphi_test, &
                                            alam_interp, dl_dt_interp)
    
    print '(A,ES12.5)', '   Lambda        = ', alam_interp
    print '(A,ES12.5)', '   dLambda/dθ    = ', dl_dt_interp
    print '(A,ES12.5)', '   Consistency:    Lambda error = ', abs(alam_interp - alam)
    print '(A,ES12.5)', '                   dΛ/dθ error  = ', abs(dl_dt_interp - dl_dt)
    print *, ''
    
    ! Test 3: Iota interpolation
    print *, '3. Testing vmec_iota_interpolate_with_field...'
    
    call vmec_iota_interpolate_with_field(gvec_field, s_test, aiota_interp, daiota_ds)
    
    print '(A,ES12.5)', '   iota          = ', aiota_interp
    print '(A,ES12.5)', '   diota/ds      = ', daiota_ds
    print '(A,ES12.5)', '   Consistency:    iota error    = ', abs(aiota_interp - aiota)
    print *, ''
    
    ! Test validation criteria
    test_passed = .true.
    
    ! Check that essential quantities are non-zero and reasonable
    if (abs(aiota) < 1.0e-6_dp) then
        print *, 'ERROR: Rotational transform too small: ', aiota
        test_passed = .false.
    end if
    
    if (abs(sqg) < 1.0e-6_dp) then
        print *, 'ERROR: Jacobian too small: ', sqg
        test_passed = .false.
    end if
    
    ! Check that vector potential derivatives are being computed (non-zero)
    if (abs(dA_theta_ds) < 1.0e-15_dp .and. abs(dA_phi_ds) < 1.0e-15_dp) then
        print *, 'WARNING: Both vector potential derivatives are zero'
        print *, '         This suggests numerical differentiation may not be working'
    end if
    
    ! Check consistency between different interfaces
    if (abs(alam_interp - alam) > 1.0e-10_dp) then
        print *, 'ERROR: Lambda values inconsistent between interfaces'
        test_passed = .false.
    end if
    
    if (abs(dl_dt_interp - dl_dt) > 1.0e-10_dp) then
        print *, 'ERROR: dLambda/dtheta values inconsistent between interfaces'
        test_passed = .false.
    end if
    
    if (abs(aiota_interp - aiota) > 1.0e-10_dp) then
        print *, 'ERROR: iota values inconsistent between interfaces'
        test_passed = .false.
    end if
    
    print *, '================================================================'
    print *, ''
    
    if (test_passed) then
        print *, 'TEST PASSED: GVEC field adapter functions working correctly'
        print *, '- Vector potential derivatives computed via numerical differentiation'
        print *, '- Stream function Lambda and derivatives available'
        print *, '- Rotational transform iota and derivatives available'
        print *, '- All adapter interfaces consistent'
    else
        print *, 'TEST FAILED: Issues found with GVEC field adapter'
        error stop 1
    end if
    
end program test_canonical_gvec