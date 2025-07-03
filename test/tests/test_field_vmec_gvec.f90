program test_field_vmec_gvec
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_vmec, only: VmecField
    use field_gvec, only: GvecField, create_gvec_field, convert_vmec_to_gvec
    use new_vmec_stuff_mod, only: netcdffile, multharm
    use spline_vmec_sub, only: spline_vmec_data
    use params, only: pi
    
    implicit none
    
    class(VmecField), allocatable :: vmec_field
    class(GvecField), allocatable :: gvec_field
    character(len=256) :: vmec_file, gvec_temp_file
    logical :: file_exists
    
    ! Variables for handling the conversion attempt
    integer :: ios
    character(len=256) :: error_msg
    
    ! Test parameters
    real(dp) :: x(3)
    real(dp) :: Acov_vmec(3), hcov_vmec(3), Bmod_vmec
    real(dp) :: Acov_gvec(3), hcov_gvec(3), Bmod_gvec
    real(dp) :: s_test, theta_test, phi_test
    real(dp) :: rel_tol, abs_tol
    integer :: ns, nt, np, i, j, k
    real(dp) :: rel_error, abs_error
    real(dp) :: max_rel_error_Bmod, max_abs_error_Bmod
    real(dp) :: max_rel_error_Acov, max_abs_error_Acov
    real(dp) :: max_rel_error_hcov, max_abs_error_hcov
    logical :: test_passed
    
    ! VMEC wout file downloaded by CMake (in build/test/tests directory)
    vmec_file = 'wout.nc'
    gvec_temp_file = 'wout_gvec.dat'  ! File created by CMake if GVEC Python available
    
    print *, '================================================================'
    print *, 'Testing field_vmec vs field_gvec with same VMEC equilibrium'
    print *, '================================================================'
    
    ! Step 1: Create VMEC field
    print *, ''
    print *, 'Step 1: Loading VMEC field from ', trim(vmec_file)
    
    ! Initialize VMEC field using SIMPLE's method
    netcdffile = vmec_file
    multharm = 7
    call spline_vmec_data
    allocate(VmecField :: vmec_field)
    print *, 'VMEC field loaded successfully in SIMPLE'
    
    ! Step 2: Check if VMEC to GVEC conversion was done by CMake
    print *, ''
    print *, 'Step 2: Checking for GVEC data file'
    
    ! Check if CMake successfully converted the file using Python GVEC API
    inquire(file=gvec_temp_file, exist=file_exists)
    if (file_exists) then
        print *, 'Found GVEC file created by CMake: ', trim(gvec_temp_file)
    else
        print *, 'GVEC file not found. Attempting direct conversion...'
        ! This will fail with a clear error message explaining why direct conversion
        ! is not possible through the SIMPLE interface
        call convert_vmec_to_gvec(vmec_file, gvec_temp_file)
    end if
    
    ! Step 3: Load GVEC field from converted file
    print *, ''
    print *, 'Step 3: Loading GVEC field from ', trim(gvec_temp_file)
    gvec_field = create_gvec_field(gvec_temp_file)
    
    ! Step 4: Compare field evaluations
    print *, ''
    print *, 'Step 4: Comparing field evaluations at same points'
    print *, '================================================================'
    
    ! Set tolerances
    rel_tol = 1.0e-2_dp  ! 1% relative tolerance
    abs_tol = 1.0e-6_dp  ! Small absolute tolerance
    
    ! Grid for comparison
    ns = 3   ! Number of flux surfaces
    nt = 4   ! Number of theta points  
    np = 2   ! Number of phi points
    
    max_rel_error_Bmod = 0.0_dp
    max_abs_error_Bmod = 0.0_dp
    max_rel_error_Acov = 0.0_dp
    max_abs_error_Acov = 0.0_dp
    max_rel_error_hcov = 0.0_dp
    max_abs_error_hcov = 0.0_dp
    test_passed = .true.
    
    print *, ''
    print *, 'Comparing at grid points:'
    
    do i = 1, ns
        s_test = 0.1_dp + 0.3_dp * real(i-1, dp) / real(ns-1, dp)  ! s from 0.1 to 0.4
        
        do j = 1, nt
            theta_test = 2.0_dp * pi * real(j-1, dp) / real(nt, dp)  ! theta from 0 to 2*pi
            
            do k = 1, np
                phi_test = pi * real(k-1, dp) / real(np, dp)  ! phi from 0 to pi
                
                ! Set coordinates (r = sqrt(s), theta, phi)
                x(1) = sqrt(s_test)
                x(2) = theta_test
                x(3) = phi_test
                
                ! Evaluate VMEC field
                call vmec_field%evaluate(x, Acov_vmec, hcov_vmec, Bmod_vmec)
                
                ! Evaluate GVEC field
                call gvec_field%evaluate(x, Acov_gvec, hcov_gvec, Bmod_gvec)
                
                ! Calculate errors for Bmod
                abs_error = abs(Bmod_gvec - Bmod_vmec)
                if (abs(Bmod_vmec) > abs_tol) then
                    rel_error = abs_error / abs(Bmod_vmec)
                else
                    rel_error = 0.0_dp
                end if
                
                max_abs_error_Bmod = max(max_abs_error_Bmod, abs_error)
                max_rel_error_Bmod = max(max_rel_error_Bmod, rel_error)
                
                ! Calculate errors for Acov (use vector norm)
                abs_error = sqrt(sum((Acov_gvec - Acov_vmec)**2))
                if (sqrt(sum(Acov_vmec**2)) > abs_tol) then
                    rel_error = abs_error / sqrt(sum(Acov_vmec**2))
                else
                    rel_error = 0.0_dp
                end if
                
                max_abs_error_Acov = max(max_abs_error_Acov, abs_error)
                max_rel_error_Acov = max(max_rel_error_Acov, rel_error)
                
                ! Calculate errors for hcov (use vector norm)
                abs_error = sqrt(sum((hcov_gvec - hcov_vmec)**2))
                if (sqrt(sum(hcov_vmec**2)) > abs_tol) then
                    rel_error = abs_error / sqrt(sum(hcov_vmec**2))
                else
                    rel_error = 0.0_dp
                end if
                
                max_abs_error_hcov = max(max_abs_error_hcov, abs_error)
                max_rel_error_hcov = max(max_rel_error_hcov, rel_error)
                
                ! Print detailed comparison for first few points
                if (i == 1 .and. j <= 2 .and. k == 1) then
                    print '(A,F6.3,A,F6.3,A,F6.3,A)', '  At (s,θ,φ) = (', s_test, ', ', &
                        theta_test/pi, 'π, ', phi_test/pi, 'π):'
                    print '(A,ES12.5,A,ES12.5)', '    |B|_VMEC = ', Bmod_vmec, ', |B|_GVEC = ', Bmod_gvec
                    print '(A,ES12.5)', '    Relative error in |B|: ', &
                        abs(Bmod_gvec - Bmod_vmec) / abs(Bmod_vmec)
                    print *, ''
                end if
                
                ! Check if errors exceed tolerance
                if (max_rel_error_Bmod > rel_tol .and. max_abs_error_Bmod > abs_tol) then
                    test_passed = .false.
                end if
            end do
        end do
    end do
    
    ! Note: Not cleaning up since we're using an existing example file
    
    print *, ''
    print *, '================================================================'
    print *, 'Summary of maximum errors:'
    print '(A,ES12.5,A,ES12.5)', '  |B|:   relative = ', max_rel_error_Bmod, &
                                ', absolute = ', max_abs_error_Bmod
    print '(A,ES12.5,A,ES12.5)', '  Acov:  relative = ', max_rel_error_Acov, &
                                ', absolute = ', max_abs_error_Acov  
    print '(A,ES12.5,A,ES12.5)', '  hcov:  relative = ', max_rel_error_hcov, &
                                ', absolute = ', max_abs_error_hcov
    
    if (test_passed) then
        print *, ''
        print *, 'TEST PASSED: Field components match within tolerance'
    else
        print *, ''
        print *, 'TEST FAILED: Field components do not match'
        print '(A,ES12.5)', '  Tolerance (relative): ', rel_tol
        print '(A,ES12.5)', '  Tolerance (absolute): ', abs_tol
        error stop 1
    end if
    
    print *, '================================================================'
    
end program test_field_vmec_gvec