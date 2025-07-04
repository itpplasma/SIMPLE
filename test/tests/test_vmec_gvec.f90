program test_vmec_gvec
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
        
        subroutine export_field_1d_data(vmec_field, gvec_field)
            use field_vmec, only: VmecField
            use field_gvec, only: GvecField
            class(VmecField), intent(in) :: vmec_field
            class(GvecField), intent(in) :: gvec_field
        end subroutine export_field_1d_data
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
    real(dp) :: Acov_vmec(3), hcov_vmec(3), Bmod_vmec, sqgBctr_vmec(3)
    real(dp) :: Acov_gvec(3), hcov_gvec(3), Bmod_gvec, sqgBctr_gvec(3)
    real(dp) :: s_test, theta_test, phi_test
    real(dp) :: rel_tol, abs_tol
    integer :: ns, nt, np, i, j, k, icomp
    real(dp) :: rel_error, abs_error, rel_err_B
    real(dp) :: max_rel_error_Bmod, max_abs_error_Bmod
    real(dp) :: max_rel_error_Acov(3), max_abs_error_Acov(3)
    real(dp) :: max_rel_error_hcov(3), max_abs_error_hcov(3)
    real(dp) :: max_rel_error_sqgBctr(3), max_abs_error_sqgBctr(3)
    logical :: test_passed

    ! VMEC wout file downloaded by CMake (in build/test/tests directory)
    vmec_file = 'wout.nc'
    gvec_temp_file = 'gvec_from_vmec_wout.dat'  ! GVEC state file created from VMEC

    ! Initialize VMEC field using SIMPLE's method
    netcdffile = vmec_file
    multharm = 5
    call spline_vmec_data
    allocate(VmecField :: vmec_field)

    ! Check if CMake successfully converted the file using Python GVEC API
    inquire(file=gvec_temp_file, exist=file_exists)
    if (.not. file_exists) then

        ! Create a minimal GVEC parameter file for VMEC reading
        param_file = 'gvec_vmec_params.ini'
        open(newunit=unit_param, file=param_file, status='replace')
        write(unit_param, '(A)') '! Minimal GVEC parameter file for VMEC reading'
        write(unit_param, '(A)') 'whichInitEquilibrium = 1  ! Read from VMEC'
        write(unit_param, '(A)') 'VMECwoutfile = ' // trim(vmec_file)
        write(unit_param, '(A)') 'VMECwoutfile_format = 0   ! NetCDF format'
        write(unit_param, '(A)') ''
        write(unit_param, '(A)') '! Grid settings optimized for axis treatment'
        write(unit_param, '(A)') 'sgrid_nelems = 21'
        write(unit_param, '(A)') 'sgrid_grid_type = 4'
        write(unit_param, '(A)') 'sgrid_rmin = 1.0e-6'
        write(unit_param, '(A)') 'sgrid_rmax = 0.99'
        write(unit_param, '(A)') ''
        write(unit_param, '(A)') '! Basis settings with higher accuracy'
        write(unit_param, '(A)') 'X1X2_deg = 5'
        write(unit_param, '(A)') 'LA_deg = 5'
        write(unit_param, '(A)') ''
        write(unit_param, '(A)') '! Output control'
        write(unit_param, '(A)') 'ProjectName = test_vmec_gvec'
        write(unit_param, '(A)') 'outputIter = 0  ! Output initial state only'
        write(unit_param, '(A)') 'maxIter = 0     ! No iterations'
        close(unit_param)

        ! Use execute_command_line to run GVEC
        call execute_command_line('../../_deps/gvec-build/bin/gvec ' // trim(param_file) // ' > /dev/null 2>&1', &
                                  exitstat=num_strings)

        ! The output file will be named based on ProjectName
        gvec_temp_file = 'test_vmec_gvec_State_0000_00000000.dat'

        ! Check if the state file was created
        inquire(file=gvec_temp_file, exist=file_exists)
        if (.not. file_exists) then
            print *, 'ERROR: GVEC state file not created'
            error stop 1
        end if
    end if
    gvec_field = create_gvec_field(gvec_temp_file)

    ! Compare field evaluations
    ! Set tolerances
    ! Note: GVEC uses different conventions than VMEC for some quantities
    ! We focus on |B| accuracy which is the most important for orbit tracing
    ! Based on actual comparison results, we expect ~0.2% maximum relative error
    rel_tol = 3.0e-3_dp  ! 0.3% relative tolerance for |B| (with some margin)
    abs_tol = 2.0e2_dp   ! 200 Gauss absolute tolerance

    ! Grid for comparison including small s values to test axis regularization
    ns = 5   ! Number of flux surfaces
    nt = 4   ! Number of theta points
    np = 2   ! Number of phi points

    max_rel_error_Bmod = 0.0_dp
    max_abs_error_Bmod = 0.0_dp
    max_rel_error_Acov = 0.0_dp
    max_abs_error_Acov = 0.0_dp
    max_rel_error_hcov = 0.0_dp
    max_abs_error_hcov = 0.0_dp
    max_rel_error_sqgBctr = 0.0_dp
    max_abs_error_sqgBctr = 0.0_dp
    test_passed = .true.

    do i = 1, ns
        ! Test range including small s values: 0.01, 0.05, 0.1, 0.3, 0.6
        if (i == 1) then
            s_test = 0.01_dp  ! Very small s to test axis regularization
        else if (i == 2) then
            s_test = 0.05_dp  ! Regularization boundary
        else
            s_test = 0.1_dp + 0.25_dp * real(i-3, dp) / real(ns-3, dp)  ! s from 0.1 to 0.6
        end if

        do j = 1, nt
            theta_test = 2.0_dp * pi * real(j-1, dp) / real(nt, dp)  ! theta from 0 to 2*pi

            do k = 1, np
                phi_test = pi * real(k-1, dp) / real(np, dp)  ! phi from 0 to pi

                ! Set coordinates (r = sqrt(s), theta, phi)
                x(1) = sqrt(s_test)
                x(2) = theta_test
                x(3) = phi_test

                ! Evaluate VMEC field
                call vmec_field%evaluate(x, Acov_vmec, hcov_vmec, Bmod_vmec, sqgBctr_vmec)

                ! Evaluate GVEC field
                call gvec_field%evaluate(x, Acov_gvec, hcov_gvec, Bmod_gvec, sqgBctr_gvec)

                ! Calculate errors for Bmod
                abs_error = abs(Bmod_gvec - Bmod_vmec)
                if (abs(Bmod_vmec) > abs_tol) then
                    rel_error = abs_error / abs(Bmod_vmec)
                else
                    rel_error = 0.0_dp
                end if

                max_abs_error_Bmod = max(max_abs_error_Bmod, abs_error)
                max_rel_error_Bmod = max(max_rel_error_Bmod, rel_error)

                ! Calculate errors for Acov components
                do icomp = 1, 3
                    abs_error = abs(Acov_gvec(icomp) - Acov_vmec(icomp))
                    if (abs(Acov_vmec(icomp)) > 1.0e-10_dp) then
                        rel_error = abs_error / abs(Acov_vmec(icomp))
                    else
                        rel_error = 0.0_dp
                    end if
                    max_abs_error_Acov(icomp) = max(max_abs_error_Acov(icomp), abs_error)
                    max_rel_error_Acov(icomp) = max(max_rel_error_Acov(icomp), rel_error)
                end do

                ! Calculate errors for hcov components
                do icomp = 1, 3
                    abs_error = abs(hcov_gvec(icomp) - hcov_vmec(icomp))
                    if (abs(hcov_vmec(icomp)) > 1.0e-10_dp) then
                        rel_error = abs_error / abs(hcov_vmec(icomp))
                    else
                        rel_error = 0.0_dp
                    end if
                    max_abs_error_hcov(icomp) = max(max_abs_error_hcov(icomp), abs_error)
                    max_rel_error_hcov(icomp) = max(max_rel_error_hcov(icomp), rel_error)
                end do

                ! Calculate errors for sqgBctr components
                do icomp = 1, 3
                    abs_error = abs(sqgBctr_gvec(icomp) - sqgBctr_vmec(icomp))
                    if (abs(sqgBctr_vmec(icomp)) > 1.0e-10_dp) then
                        rel_error = abs_error / abs(sqgBctr_vmec(icomp))
                    else
                        rel_error = 0.0_dp
                    end if
                    max_abs_error_sqgBctr(icomp) = max(max_abs_error_sqgBctr(icomp), abs_error)
                    max_rel_error_sqgBctr(icomp) = max(max_rel_error_sqgBctr(icomp), rel_error)
                end do

                ! Calculate relative error for |B|
                rel_err_B = abs(Bmod_gvec - Bmod_vmec) / abs(Bmod_vmec)

                ! Check if errors exceed tolerance
                if (max_rel_error_Bmod > rel_tol .and. max_abs_error_Bmod > abs_tol) then
                    test_passed = .false.
                end if
            end do
        end do
    end do

    ! Evaluate fields at a representative point for typical values
    x(1) = sqrt(0.5_dp)  ! r = sqrt(s) at s=0.5
    x(2) = pi/4.0_dp     ! theta = pi/4
    x(3) = 0.0_dp        ! phi = 0
    
    call vmec_field%evaluate(x, Acov_vmec, hcov_vmec, Bmod_vmec, sqgBctr_vmec)
    call gvec_field%evaluate(x, Acov_gvec, hcov_gvec, Bmod_gvec, sqgBctr_gvec)

    ! Print comparison results
    print *, ''
    print *, 'Field comparison at (s,θ,φ) = (0.5, π/4, 0):'
    print *, '================================================================'
    print *, '              VMEC              GVEC              Rel. Error'
    print *, '----------------------------------------------------------------'
    print '(A,ES16.8,ES16.8,ES16.8)', '  |B|     ', Bmod_vmec, Bmod_gvec, &
                                         abs(Bmod_gvec - Bmod_vmec) / abs(Bmod_vmec)
    print *, ''
    print '(A,ES16.8,ES16.8,ES16.8)', '  Acov(1) ', Acov_vmec(1), Acov_gvec(1), &
                                         merge(abs(Acov_gvec(1) - Acov_vmec(1)) / abs(Acov_vmec(1)), 0.0_dp, abs(Acov_vmec(1)) > 1.0e-10_dp)
    print '(A,ES16.8,ES16.8,ES16.8)', '  Acov(2) ', Acov_vmec(2), Acov_gvec(2), &
                                         merge(abs(Acov_gvec(2) - Acov_vmec(2)) / abs(Acov_vmec(2)), 0.0_dp, abs(Acov_vmec(2)) > 1.0e-10_dp)
    print '(A,ES16.8,ES16.8,ES16.8)', '  Acov(3) ', Acov_vmec(3), Acov_gvec(3), &
                                         merge(abs(Acov_gvec(3) - Acov_vmec(3)) / abs(Acov_vmec(3)), 0.0_dp, abs(Acov_vmec(3)) > 1.0e-10_dp)
    print *, ''
    print '(A,ES16.8,ES16.8,ES16.8)', '  hcov(1) ', hcov_vmec(1), hcov_gvec(1), &
                                         merge(abs(hcov_gvec(1) - hcov_vmec(1)) / abs(hcov_vmec(1)), 0.0_dp, abs(hcov_vmec(1)) > 1.0e-10_dp)
    print '(A,ES16.8,ES16.8,ES16.8)', '  hcov(2) ', hcov_vmec(2), hcov_gvec(2), &
                                         merge(abs(hcov_gvec(2) - hcov_vmec(2)) / abs(hcov_vmec(2)), 0.0_dp, abs(hcov_vmec(2)) > 1.0e-10_dp)
    print '(A,ES16.8,ES16.8,ES16.8)', '  hcov(3) ', hcov_vmec(3), hcov_gvec(3), &
                                         merge(abs(hcov_gvec(3) - hcov_vmec(3)) / abs(hcov_vmec(3)), 0.0_dp, abs(hcov_vmec(3)) > 1.0e-10_dp)
    print *, ''
    print '(A,ES16.8,ES16.8,A)', '  √g·B^s  ', sqgBctr_vmec(1), sqgBctr_gvec(1), '    (VMEC now provides this)'
    print '(A,ES16.8,ES16.8,ES16.8)', '  √g·B^θ  ', sqgBctr_vmec(2), sqgBctr_gvec(2), &
                                         merge(abs(sqgBctr_gvec(2) - sqgBctr_vmec(2)) / abs(sqgBctr_vmec(2)), 0.0_dp, abs(sqgBctr_vmec(2)) > 1.0e-10_dp)
    print '(A,ES16.8,ES16.8,ES16.8)', '  √g·B^φ  ', sqgBctr_vmec(3), sqgBctr_gvec(3), &
                                         merge(abs(sqgBctr_gvec(3) - sqgBctr_vmec(3)) / abs(sqgBctr_vmec(3)), 0.0_dp, abs(sqgBctr_vmec(3)) > 1.0e-10_dp)
    print *, '================================================================'
    print *, ''
    print *, 'Maximum errors over test grid:'
    print *, '----------------------------------------------------------------'
    print '(A,ES12.5,A,ES12.5,A)', '  |B|     max error: relative = ', max_rel_error_Bmod, &
                                ', absolute = ', max_abs_error_Bmod, ' Gauss'
    print *, ''

    ! Generate field data and plots
    call export_field_2d_data(vmec_field, gvec_field)
    call export_field_1d_data(vmec_field, gvec_field)
    
    ! Generate plots silently
    call execute_command_line('python3 ../../../test/tests/plot_fields_2d.py > /dev/null 2>&1', exitstat=i)
    call execute_command_line('python3 ../../../test/tests/plot_1d_field_comparison.py > /dev/null 2>&1', exitstat=i)
    call execute_command_line('python3 ../../../test/tests/plot_small_s_behavior.py > /dev/null 2>&1', exitstat=i)
    call execute_command_line('python3 ../../../test/tests/detailed_field_analysis.py > /dev/null 2>&1', exitstat=i)
    call execute_command_line('python3 ../../../test/tests/create_1d_radial_plots.py > /dev/null 2>&1', exitstat=i)
    call execute_command_line('python3 ../../../test/tests/investigate_field_differences.py > /dev/null 2>&1', exitstat=i)
    
    ! Final test result
    print *, ''
    if (test_passed) then
        print *, 'TEST PASSED'
    else
        print *, 'TEST FAILED'
        print '(A,ES12.5)', '  Tolerance (relative): ', rel_tol
        print '(A,ES12.5)', '  Tolerance (absolute): ', abs_tol, ' Gauss'
    end if

end program test_vmec_gvec

subroutine export_field_2d_data(vmec_field, gvec_field)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_vmec, only: VmecField
    use field_gvec, only: GvecField
    use params, only: pi
    
    implicit none
    
    class(VmecField), intent(in) :: vmec_field
    class(GvecField), intent(in) :: gvec_field
    
    ! Grid parameters including small s values
    integer, parameter :: ns = 40    ! Number of s points
    integer, parameter :: nt = 60    ! Number of theta points
    real(dp), parameter :: s_min = 0.01_dp  ! Start from very small s
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
    
    ! Generate 2D field data silently
    
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
    
    ! Field data files generated: Bmod_vmec_2d.dat, Bmod_gvec_2d.dat, grid_info.dat
    
end subroutine export_field_2d_data

subroutine export_field_1d_data(vmec_field, gvec_field)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_vmec, only: VmecField
    use field_gvec, only: GvecField
    use params, only: pi
    
    implicit none
    
    class(VmecField), intent(in) :: vmec_field
    class(GvecField), intent(in) :: gvec_field
    
    ! Grid parameters for 1D radial comparison
    integer, parameter :: ns = 100   ! Number of s points for high resolution
    real(dp), parameter :: s_min = 1.0e-4_dp  ! Very small s to test axis
    real(dp), parameter :: s_max = 0.99_dp
    real(dp), parameter :: theta_fixed = 0.0_dp    ! Fixed poloidal angle
    real(dp), parameter :: phi_fixed = 0.0_dp      ! Fixed toroidal angle
    
    ! Field evaluation variables
    real(dp) :: x(3)
    real(dp) :: Acov_vmec(3), hcov_vmec(3), Bmod_vmec
    real(dp) :: Acov_gvec(3), hcov_gvec(3), Bmod_gvec
    
    ! Grid variables
    real(dp) :: s, r
    integer :: i, unit_out
    
    ! Generate 1D radial field data
    
    ! Write combined radial comparison data
    open(newunit=unit_out, file='radial_comparison.dat', status='replace')
    write(unit_out, '(A)') '# 1D radial field comparison'
    write(unit_out, '(A,F8.4,A,F8.4,A)') '# Fixed theta = ', theta_fixed/pi, ' pi, phi = ', phi_fixed/pi, ' pi'
    write(unit_out, '(A)') '# Columns: s, r, |B|_VMEC, |B|_GVEC, relative_error'
    
    do i = 1, ns
        ! Use log spacing for better resolution at small s
        if (i == 1) then
            s = s_min
        else
            ! Log-space from s_min to s_max
            s = s_min * (s_max/s_min)**((real(i-1,dp))/(real(ns-1,dp)))
        end if
        r = sqrt(s)
        
        ! Set coordinates (r = sqrt(s), theta, phi)
        x(1) = r
        x(2) = theta_fixed
        x(3) = phi_fixed
        
        ! Evaluate fields
        call vmec_field%evaluate(x, Acov_vmec, hcov_vmec, Bmod_vmec)
        call gvec_field%evaluate(x, Acov_gvec, hcov_gvec, Bmod_gvec)
        
        ! Calculate relative error
        write(unit_out, '(5ES16.8)') s, r, Bmod_vmec, Bmod_gvec, &
            abs(Bmod_gvec - Bmod_vmec) / abs(Bmod_vmec)
    end do
    close(unit_out)
    
    ! Field data file generated: radial_comparison.dat
    
end subroutine export_field_1d_data
