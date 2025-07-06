program export_field_2d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_vmec, only: VmecField
    use field_gvec, only: GvecField, create_gvec_field
    use new_vmec_stuff_mod, only: netcdffile, multharm
    use spline_vmec_sub, only: spline_vmec_data
    use params, only: pi
    
    implicit none
    
    class(VmecField), allocatable :: vmec_field
    class(GvecField), allocatable :: gvec_field
    character(len=256) :: vmec_file, gvec_file
    
    ! Grid parameters
    integer, parameter :: ns = 50    ! Number of s points
    integer, parameter :: nt = 100   ! Number of theta points
    real(dp) :: s_min = 0.1_dp
    real(dp) :: s_max = 0.9_dp
    real(dp) :: theta_min = 0.0_dp
    real(dp) :: theta_max = 2.0_dp * pi
    real(dp) :: phi_fixed = 0.0_dp  ! Fixed toroidal angle for 2D slice
    
    ! Field evaluation variables
    real(dp) :: x(3)
    real(dp) :: Acov_vmec(3), hcov_vmec(3), Bmod_vmec
    real(dp) :: Acov_gvec(3), hcov_gvec(3), Bmod_gvec
    
    ! Grid variables
    real(dp), allocatable :: s_grid(:), theta_grid(:)
    real(dp), allocatable :: Bmod_vmec_2d(:,:), Bmod_gvec_2d(:,:)
    real(dp), allocatable :: R_grid(:,:), Z_grid(:,:)
    
    integer :: i, j, unit_out
    real(dp) :: s, theta, r
    
    ! File names
    vmec_file = 'wout.nc'
    gvec_file = 'test_vmec_gvec_State_0000_00000000.dat'
    
    print *, '================================================================'
    print *, 'Exporting 2D magnetic field data for VMEC and GVEC'
    print *, '================================================================'
    
    ! Allocate grids
    allocate(s_grid(ns), theta_grid(nt))
    allocate(Bmod_vmec_2d(ns, nt), Bmod_gvec_2d(ns, nt))
    allocate(R_grid(ns, nt), Z_grid(ns, nt))
    
    ! Create grids
    do i = 1, ns
        s_grid(i) = s_min + (s_max - s_min) * real(i-1, dp) / real(ns-1, dp)
    end do
    
    do j = 1, nt
        theta_grid(j) = theta_min + (theta_max - theta_min) * real(j-1, dp) / real(nt-1, dp)
    end do
    
    ! Initialize VMEC field
    print *, ''
    print *, 'Loading VMEC field from ', trim(vmec_file)
    netcdffile = vmec_file
    multharm = 7
    call spline_vmec_data
    allocate(VmecField :: vmec_field)
    
    ! Initialize GVEC field
    print *, 'Loading GVEC field from ', trim(gvec_file)
    gvec_field = create_gvec_field(gvec_file)
    
    ! Evaluate fields on 2D grid
    print *, ''
    print *, 'Evaluating fields on 2D grid (s, theta) at phi = ', phi_fixed/pi, 'π'
    print *, 'Grid size: ', ns, ' x ', nt
    
    do i = 1, ns
        s = s_grid(i)
        r = sqrt(s)  ! r = sqrt(s)
        
        do j = 1, nt
            theta = theta_grid(j)
            
            ! Set coordinates (r = sqrt(s), theta, phi)
            x(1) = r
            x(2) = theta
            x(3) = phi_fixed
            
            ! Evaluate VMEC field
            call vmec_field%evaluate(x, Acov_vmec, hcov_vmec, Bmod_vmec)
            Bmod_vmec_2d(i, j) = Bmod_vmec
            
            ! Evaluate GVEC field
            call gvec_field%evaluate(x, Acov_gvec, hcov_gvec, Bmod_gvec)
            Bmod_gvec_2d(i, j) = Bmod_gvec
            
            ! For visualization, we'll also need R,Z coordinates
            ! Approximate R,Z for plotting (assuming axisymmetric at phi=0)
            ! This is a rough approximation - in reality we'd compute from VMEC harmonics
            R_grid(i, j) = 5.0_dp + r * cos(theta)  ! Major radius ~ 5m
            Z_grid(i, j) = r * sin(theta)
        end do
        
        if (mod(i, 10) == 0) then
            print '(A,I3,A)', '  Progress: ', nint(100.0_dp * real(i, dp) / real(ns, dp)), '%'
        end if
    end do
    
    ! Write VMEC data
    print *, ''
    print *, 'Writing VMEC field data to Bmod_vmec_2d.dat'
    open(newunit=unit_out, file='Bmod_vmec_2d.dat', status='replace')
    write(unit_out, '(A)') '# 2D magnetic field magnitude from VMEC'
    write(unit_out, '(A,F8.4,A)') '# Fixed phi = ', phi_fixed/pi, ' pi'
    write(unit_out, '(A,I5,A,I5)') '# Grid size: ns = ', ns, ', nt = ', nt
    write(unit_out, '(A)') '# Columns: s, theta, R, Z, |B|'
    
    do i = 1, ns
        do j = 1, nt
            write(unit_out, '(5ES16.8)') s_grid(i), theta_grid(j), &
                R_grid(i, j), Z_grid(i, j), Bmod_vmec_2d(i, j)
        end do
        write(unit_out, *)  ! Empty line for gnuplot
    end do
    close(unit_out)
    
    ! Write GVEC data
    print *, 'Writing GVEC field data to Bmod_gvec_2d.dat'
    open(newunit=unit_out, file='Bmod_gvec_2d.dat', status='replace')
    write(unit_out, '(A)') '# 2D magnetic field magnitude from GVEC'
    write(unit_out, '(A,F8.4,A)') '# Fixed phi = ', phi_fixed/pi, ' pi'
    write(unit_out, '(A,I5,A,I5)') '# Grid size: ns = ', ns, ', nt = ', nt
    write(unit_out, '(A)') '# Columns: s, theta, R, Z, |B|'
    
    do i = 1, ns
        do j = 1, nt
            write(unit_out, '(5ES16.8)') s_grid(i), theta_grid(j), &
                R_grid(i, j), Z_grid(i, j), Bmod_gvec_2d(i, j)
        end do
        write(unit_out, *)  ! Empty line for gnuplot
    end do
    close(unit_out)
    
    ! Write grid info for Python script
    print *, 'Writing grid information to grid_info.dat'
    open(newunit=unit_out, file='grid_info.dat', status='replace')
    write(unit_out, '(I5,I5)') ns, nt
    write(unit_out, '(4ES16.8)') s_min, s_max, theta_min, theta_max
    write(unit_out, '(ES16.8)') phi_fixed
    close(unit_out)
    
    ! Print summary statistics
    print *, ''
    print *, 'Field statistics:'
    print '(A,ES12.5,A,ES12.5)', '  VMEC |B| range: ', &
        minval(Bmod_vmec_2d), ' to ', maxval(Bmod_vmec_2d)
    print '(A,ES12.5,A,ES12.5)', '  GVEC |B| range: ', &
        minval(Bmod_gvec_2d), ' to ', maxval(Bmod_gvec_2d)
    
    ! Deallocate
    deallocate(s_grid, theta_grid)
    deallocate(Bmod_vmec_2d, Bmod_gvec_2d)
    deallocate(R_grid, Z_grid)
    
    print *, ''
    print *, 'Export completed successfully!'
    print *, 'Use plot_fields_2d.py to visualize the results'
    print *, '================================================================'
    
end program export_field_2d