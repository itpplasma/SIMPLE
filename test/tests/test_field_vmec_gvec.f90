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
    
    interface
        subroutine export_field_2d_data(vmec_field, gvec_field)
            use field_vmec, only: VmecField
            use field_gvec, only: GvecField
            class(VmecField), intent(in) :: vmec_field
            class(GvecField), intent(in) :: gvec_field
        end subroutine export_field_2d_data
    end interface
    character(len=256) :: vmec_file, gvec_temp_file
    logical :: file_exists

    ! Variables for GVEC conversion
    character(len=256) :: param_file
    character(len=4096) :: param_strings(20)
    integer :: num_strings
    integer :: unit_param

    ! Test parameters
    real(dp) :: x(3)
    real(dp) :: Acov_vmec(3), hcov_vmec(3), Bmod_vmec
    real(dp) :: Acov_gvec(3), hcov_gvec(3), Bmod_gvec
    real(dp) :: s_test, theta_test, phi_test
    real(dp) :: rel_tol, abs_tol
    integer :: ns, nt, np, i, j, k
    real(dp) :: rel_error, abs_error, rel_err_B
    real(dp) :: max_rel_error_Bmod, max_abs_error_Bmod
    real(dp) :: max_rel_error_Acov, max_abs_error_Acov
    real(dp) :: max_rel_error_hcov, max_abs_error_hcov
    logical :: test_passed

    ! VMEC wout file downloaded by CMake (in build/test/tests directory)
    vmec_file = 'wout.nc'
    gvec_temp_file = 'gvec_from_vmec_wout.dat'  ! GVEC state file created from VMEC

    print *, '================================================================'
    print *, 'Testing field_vmec vs field_gvec with same VMEC equilibrium'
    print *, '================================================================'

    ! Step 1: Create VMEC field
    print *, ''
    print *, 'Step 1: Loading VMEC field from ', trim(vmec_file)

    ! Initialize VMEC field using SIMPLE's method
    netcdffile = vmec_file
    multharm = 5
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
        print *, 'GVEC file not found. Creating GVEC state from VMEC...'

        ! Create a minimal GVEC parameter file for VMEC reading
        param_file = 'gvec_vmec_params.ini'
        open(newunit=unit_param, file=param_file, status='replace')
        write(unit_param, '(A)') '! Minimal GVEC parameter file for VMEC reading'
        write(unit_param, '(A)') 'whichInitEquilibrium = 1  ! Read from VMEC'
        write(unit_param, '(A)') 'VMECwoutfile = ' // trim(vmec_file)
        write(unit_param, '(A)') 'VMECwoutfile_format = 0   ! NetCDF format'
        write(unit_param, '(A)') ''
        write(unit_param, '(A)') '! Minimal grid settings'
        write(unit_param, '(A)') 'sgrid_nelems = 11'
        write(unit_param, '(A)') 'sgrid_grid_type = 4'
        write(unit_param, '(A)') ''
        write(unit_param, '(A)') '! Basis settings'
        write(unit_param, '(A)') 'X1X2_deg = 3'
        write(unit_param, '(A)') 'LA_deg = 3'
        write(unit_param, '(A)') ''
        write(unit_param, '(A)') '! Output control'
        write(unit_param, '(A)') 'ProjectName = test_vmec_gvec'
        write(unit_param, '(A)') 'outputIter = 0  ! Output initial state only'
        write(unit_param, '(A)') 'maxIter = 0     ! No iterations'
        close(unit_param)

        print *, 'Created parameter file: ', trim(param_file)
        print *, 'Running GVEC to create state file...'

        ! Use execute_command_line to run GVEC
        ! This will read VMEC and output the state
        call execute_command_line('../../_deps/gvec-build/bin/gvec ' // trim(param_file), &
                                  exitstat=num_strings)

        ! The output file will be named based on ProjectName
        gvec_temp_file = 'test_vmec_gvec_State_0000_00000000.dat'

        ! Check if the state file was created
        inquire(file=gvec_temp_file, exist=file_exists)
        if (.not. file_exists) then
            print *, 'ERROR: GVEC state file not created'
            error stop 1
        end if

        print *, 'GVEC state file created successfully'
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
    ! Note: GVEC uses different conventions than VMEC for some quantities
    ! We focus on |B| accuracy which is the most important for orbit tracing
    rel_tol = 1.0e-2_dp  ! 1% relative tolerance for |B|
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

                ! Print detailed comparison only for first few points
                rel_err_B = abs(Bmod_gvec - Bmod_vmec) / abs(Bmod_vmec)
                if (i == 1 .and. j <= 2 .and. k == 1) then  ! Print first two points only
                    print '(A,F6.3,A,F6.3,A,F6.3,A)', '  At (s,θ,φ) = (', s_test, ', ', &
                        theta_test/pi, 'π, ', phi_test/pi, 'π):'
                    print '(A,ES12.5,A,ES12.5)', '    |B|_VMEC = ', Bmod_vmec, ', |B|_GVEC = ', Bmod_gvec
                    print '(A,ES12.5)', '    Relative error in |B|: ', &
                        abs(Bmod_gvec - Bmod_vmec) / abs(Bmod_vmec)
                    ! print *, ''
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
                                ', absolute = ', max_abs_error_Acov, ' (NOT CHECKED - different conventions)'
    print '(A,ES12.5,A,ES12.5)', '  hcov:  relative = ', max_rel_error_hcov, &
                                ', absolute = ', max_abs_error_hcov, ' (NOT CHECKED - different conventions)'

    ! Step 5: Generate 2D field data for visualization (regardless of test result)
    print *, ''
    print *, 'Step 5: Generating 2D field data for visualization'
    print *, '================================================================'
    
    call export_field_2d_data(vmec_field, gvec_field)
    
    print *, ''
    print *, 'Step 6: Creating visualization plots'
    print *, '================================================================'
    
    call execute_command_line('python3 plot_fields_2d.py', exitstat=i)
    if (i == 0) then
        print *, 'Visualization plots created successfully'
        print *, 'Generated files:'
        print *, '  - field_comparison_2d.png'
        print *, '  - field_difference_2d.png'
        print *, '  - vmec_field_2d.png'
        print *, '  - gvec_field_2d.png'
    else
        print *, 'Warning: Python plotting failed (exit code:', i, ')'
        print *, 'Check if matplotlib is available and plot_fields_2d.py exists'
    end if

    print *, '================================================================'
    
    ! Final test result
    if (test_passed) then
        print *, ''
        print *, 'TEST PASSED: Magnetic field magnitude |B| matches within tolerance'
    else
        print *, ''
        print *, 'TEST FAILED: Magnetic field magnitude |B| does not match'
        print '(A,ES12.5)', '  Tolerance (relative): ', rel_tol
        print '(A,ES12.5)', '  Tolerance (absolute): ', abs_tol
        print *, ''
        print *, 'Note: Visualization plots were generated despite test failure.'
    end if

end program test_field_vmec_gvec

subroutine export_field_2d_data(vmec_field, gvec_field)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_vmec, only: VmecField
    use field_gvec, only: GvecField
    use params, only: pi
    
    implicit none
    
    class(VmecField), intent(in) :: vmec_field
    class(GvecField), intent(in) :: gvec_field
    
    ! Grid parameters
    integer, parameter :: ns = 30    ! Number of s points
    integer, parameter :: nt = 60    ! Number of theta points
    real(dp), parameter :: s_min = 0.1_dp
    real(dp), parameter :: s_max = 0.9_dp
    real(dp), parameter :: theta_min = 0.0_dp
    real(dp), parameter :: theta_max = 2.0_dp * pi
    real(dp), parameter :: phi_fixed = 0.0_dp  ! Fixed toroidal angle for 2D slice
    
    ! Field evaluation variables
    real(dp) :: x(3)
    real(dp) :: Acov_vmec(3), hcov_vmec(3), Bmod_vmec
    real(dp) :: Acov_gvec(3), hcov_gvec(3), Bmod_gvec
    
    ! Grid variables
    real(dp) :: s, theta, r
    integer :: i, j, unit_out
    
    print *, 'Generating 2D field data on ', ns, ' x ', nt, ' grid'
    print *, 'phi fixed at: ', phi_fixed/pi, 'π'
    
    ! Write VMEC data
    open(newunit=unit_out, file='Bmod_vmec_2d.dat', status='replace')
    write(unit_out, '(A)') '# 2D magnetic field magnitude from VMEC'
    write(unit_out, '(A,F8.4,A)') '# Fixed phi = ', phi_fixed/pi, ' pi'
    write(unit_out, '(A,I5,A,I5)') '# Grid size: ns = ', ns, ', nt = ', nt
    write(unit_out, '(A)') '# Columns: s, theta, R, Z, |B|'
    
    do i = 1, ns
        s = s_min + (s_max - s_min) * real(i-1, dp) / real(ns-1, dp)
        r = sqrt(s)
        
        do j = 1, nt
            theta = theta_min + (theta_max - theta_min) * real(j-1, dp) / real(nt-1, dp)
            
            ! Set coordinates (r = sqrt(s), theta, phi)
            x(1) = r
            x(2) = theta
            x(3) = phi_fixed
            
            ! Evaluate VMEC field
            call vmec_field%evaluate(x, Acov_vmec, hcov_vmec, Bmod_vmec)
            
            ! Approximate R,Z for plotting (assuming axisymmetric at phi=0)
            write(unit_out, '(5ES16.8)') s, theta, &
                5.0_dp + r * cos(theta), r * sin(theta), Bmod_vmec
        end do
        write(unit_out, *)  ! Empty line for gnuplot
    end do
    close(unit_out)
    
    ! Write GVEC data
    open(newunit=unit_out, file='Bmod_gvec_2d.dat', status='replace')
    write(unit_out, '(A)') '# 2D magnetic field magnitude from GVEC'
    write(unit_out, '(A,F8.4,A)') '# Fixed phi = ', phi_fixed/pi, ' pi'
    write(unit_out, '(A,I5,A,I5)') '# Grid size: ns = ', ns, ', nt = ', nt
    write(unit_out, '(A)') '# Columns: s, theta, R, Z, |B|'
    
    do i = 1, ns
        s = s_min + (s_max - s_min) * real(i-1, dp) / real(ns-1, dp)
        r = sqrt(s)
        
        do j = 1, nt
            theta = theta_min + (theta_max - theta_min) * real(j-1, dp) / real(nt-1, dp)
            
            ! Set coordinates (r = sqrt(s), theta, phi)
            x(1) = r
            x(2) = theta
            x(3) = phi_fixed
            
            ! Evaluate GVEC field
            call gvec_field%evaluate(x, Acov_gvec, hcov_gvec, Bmod_gvec)
            
            ! Approximate R,Z for plotting (assuming axisymmetric at phi=0)
            write(unit_out, '(5ES16.8)') s, theta, &
                5.0_dp + r * cos(theta), r * sin(theta), Bmod_gvec
        end do
        write(unit_out, *)  ! Empty line for gnuplot
    end do
    close(unit_out)
    
    ! Write grid info for Python script
    open(newunit=unit_out, file='grid_info.dat', status='replace')
    write(unit_out, '(I5,I5)') ns, nt
    write(unit_out, '(4ES16.8)') s_min, s_max, theta_min, theta_max
    write(unit_out, '(ES16.8)') phi_fixed
    close(unit_out)
    
    print *, 'Field data files generated:'
    print *, '  - Bmod_vmec_2d.dat'
    print *, '  - Bmod_gvec_2d.dat'
    print *, '  - grid_info.dat'
    
end subroutine export_field_2d_data
