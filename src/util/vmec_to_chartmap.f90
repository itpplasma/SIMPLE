program vmec_to_chartmap
    !> Generate chartmap coordinate file from VMEC equilibrium.
    !>
    !> Two modes:
    !>   --passthrough: Use libneo VMEC splines directly (exact match)
    !>   --map2disc: Call Python map2disc for conformal mapping
    !>
    !> Usage:
    !>   vmec_to_chartmap.x wout.nc output.chartmap.nc --passthrough [--nrho N] [--ntheta N] [--nzeta N]
    !>   vmec_to_chartmap.x wout.nc output.chartmap.nc --map2disc [--nrho N] [--ntheta N] [--nzeta N]

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use libneo_coordinates, only: coordinate_system_t, make_vmec_coordinate_system
    use util, only: twopi
    use netcdf

    implicit none

    character(256) :: vmec_file, output_file, mode, arg
    integer :: nrho, ntheta, nzeta, nfp
    integer :: i, nargs
    real(dp) :: dummy
    logical :: use_map2disc

    ! Defaults
    nrho = 64
    ntheta = 65
    nzeta = 66
    use_map2disc = .false.

    nargs = command_argument_count()
    if (nargs < 3) then
        call print_usage()
        stop 1
    end if

    call get_command_argument(1, vmec_file)
    call get_command_argument(2, output_file)
    call get_command_argument(3, mode)

    if (trim(mode) == '--passthrough') then
        use_map2disc = .false.
    else if (trim(mode) == '--map2disc') then
        use_map2disc = .true.
    else
        print *, 'Error: mode must be --passthrough or --map2disc'
        call print_usage()
        stop 1
    end if

    ! Parse optional arguments
    i = 4
    do while (i <= nargs)
        call get_command_argument(i, arg)
        if (trim(arg) == '--nrho') then
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) nrho
        else if (trim(arg) == '--ntheta') then
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) ntheta
        else if (trim(arg) == '--nzeta') then
            i = i + 1
            call get_command_argument(i, arg)
            read(arg, *) nzeta
        else
            print *, 'Unknown argument: ', trim(arg)
            call print_usage()
            stop 1
        end if
        i = i + 1
    end do

    if (use_map2disc) then
        call generate_map2disc(vmec_file, output_file, nrho, ntheta, nzeta)
    else
        call generate_passthrough(vmec_file, output_file, nrho, ntheta, nzeta)
    end if

contains

    subroutine print_usage()
        print *, 'Usage: vmec_to_chartmap.x <vmec_file> <output_file> <mode> [options]'
        print *, ''
        print *, 'Modes:'
        print *, '  --passthrough  Use libneo VMEC splines directly'
        print *, '  --map2disc     Use Python map2disc conformal mapping'
        print *, ''
        print *, 'Options:'
        print *, '  --nrho N       Number of radial points (default: 64)'
        print *, '  --ntheta N     Number of poloidal points (default: 65)'
        print *, '  --nzeta N      Number of toroidal points (default: 66)'
    end subroutine print_usage


    subroutine generate_passthrough(vmec_file, output_file, nrho, ntheta, nzeta)
        !> Generate chartmap using libneo VMEC splines directly.
        character(*), intent(in) :: vmec_file, output_file
        integer, intent(in) :: nrho, ntheta, nzeta

        class(coordinate_system_t), allocatable :: vmec_coords
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
        real(dp) :: u(3), x_cyl(3)
        integer :: ir, ith, iz, nfp
        real(dp) :: dummy

        print *, 'Generating chartmap using libneo VMEC splines (passthrough)'
        print *, '  VMEC file: ', trim(vmec_file)
        print *, '  Output file: ', trim(output_file)
        print *, '  Grid: ', nrho, ' x ', ntheta, ' x ', nzeta

        ! Initialize VMEC
        call init_vmec(vmec_file, 5, 5, 5, dummy)
        call make_vmec_coordinate_system(vmec_coords)

        ! Get nfp from VMEC file
        call get_vmec_nfp(vmec_file, nfp)

        ! Create coordinate grids
        allocate(rho(nrho), theta(ntheta), zeta(nzeta))
        allocate(x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, ntheta, nzeta))

        do ir = 1, nrho
            rho(ir) = real(ir - 1, dp) / real(nrho - 1, dp)
        end do
        do ith = 1, ntheta
            theta(ith) = twopi * real(ith - 1, dp) / real(ntheta, dp)
        end do
        do iz = 1, nzeta
            zeta(iz) = (twopi / nfp) * real(iz - 1, dp) / real(nzeta, dp)
        end do

        ! Evaluate VMEC coordinates and convert to Cartesian
        do iz = 1, nzeta
            do ith = 1, ntheta
                do ir = 1, nrho
                    u(1) = rho(ir)
                    u(2) = theta(ith)
                    u(3) = zeta(iz)

                    ! VMEC returns cylindrical (R, phi, Z)
                    call vmec_coords%evaluate_point(u, x_cyl)

                    ! Convert to Cartesian
                    x(ir, ith, iz) = x_cyl(1) * cos(x_cyl(2))
                    y(ir, ith, iz) = x_cyl(1) * sin(x_cyl(2))
                    z(ir, ith, iz) = x_cyl(3)
                end do
            end do
        end do

        ! Write chartmap file
        call write_chartmap_nc(output_file, rho, theta, zeta, x, y, z, nfp)

        deallocate(vmec_coords)
        deallocate(rho, theta, zeta, x, y, z)

        print *, 'Done.'
    end subroutine generate_passthrough


    subroutine generate_map2disc(vmec_file, output_file, nrho, ntheta, nzeta)
        !> Generate chartmap using Python map2disc.
        character(*), intent(in) :: vmec_file, output_file
        integer, intent(in) :: nrho, ntheta, nzeta

        character(1024) :: cmd
        integer :: exitstat, cmdstat

        print *, 'Generating chartmap using Python map2disc'
        print *, '  VMEC file: ', trim(vmec_file)
        print *, '  Output file: ', trim(output_file)
        print *, '  Grid: ', nrho, ' x ', ntheta, ' x ', nzeta

        write(cmd, '(A,A,A,A,A,I0,A,I0,A,I0)') &
            'python3 -c "from simple.chartmap import vmec_to_chartmap_map2disc; ', &
            'vmec_to_chartmap_map2disc(''', trim(vmec_file), ''', ''', trim(output_file), &
            ''', nrho=', nrho, ', ntheta=', ntheta, ', nzeta=', nzeta, ')"'

        print *, '  Command: ', trim(cmd)

        call execute_command_line(trim(cmd), wait=.true., exitstat=exitstat, cmdstat=cmdstat)

        if (cmdstat /= 0 .or. exitstat /= 0) then
            print *, 'Error: map2disc failed with exit status ', exitstat
            stop 1
        end if

        print *, 'Done.'
    end subroutine generate_map2disc


    subroutine get_vmec_nfp(filename, nfp)
        character(*), intent(in) :: filename
        integer, intent(out) :: nfp

        integer :: ncid, varid, ierr

        ierr = nf90_open(filename, NF90_NOWRITE, ncid)
        if (ierr /= NF90_NOERR) then
            print *, 'Error opening VMEC file: ', trim(filename)
            stop 1
        end if

        ierr = nf90_inq_varid(ncid, 'nfp', varid)
        ierr = nf90_get_var(ncid, varid, nfp)
        ierr = nf90_close(ncid)
    end subroutine get_vmec_nfp


    subroutine write_chartmap_nc(filename, rho, theta, zeta, x, y, z, nfp)
        character(*), intent(in) :: filename
        real(dp), intent(in) :: rho(:), theta(:), zeta(:)
        real(dp), intent(in) :: x(:,:,:), y(:,:,:), z(:,:,:)
        integer, intent(in) :: nfp

        integer :: ncid, rho_dimid, theta_dimid, zeta_dimid
        integer :: rho_varid, theta_varid, zeta_varid
        integer :: x_varid, y_varid, z_varid, nfp_varid
        integer :: ierr, nrho, ntheta, nzeta

        nrho = size(rho)
        ntheta = size(theta)
        nzeta = size(zeta)

        ! Create file
        ierr = nf90_create(filename, NF90_CLOBBER, ncid)
        if (ierr /= NF90_NOERR) then
            print *, 'Error creating file: ', trim(filename)
            stop 1
        end if

        ! Define dimensions
        ierr = nf90_def_dim(ncid, 'rho', nrho, rho_dimid)
        ierr = nf90_def_dim(ncid, 'theta', ntheta, theta_dimid)
        ierr = nf90_def_dim(ncid, 'zeta', nzeta, zeta_dimid)

        ! Define variables
        ierr = nf90_def_var(ncid, 'rho', NF90_DOUBLE, (/rho_dimid/), rho_varid)
        ierr = nf90_def_var(ncid, 'theta', NF90_DOUBLE, (/theta_dimid/), theta_varid)
        ierr = nf90_def_var(ncid, 'zeta', NF90_DOUBLE, (/zeta_dimid/), zeta_varid)
        ! Note: libneo expects (zeta, theta, rho) ordering in file
        ierr = nf90_def_var(ncid, 'x', NF90_DOUBLE, &
            (/rho_dimid, theta_dimid, zeta_dimid/), x_varid)
        ierr = nf90_def_var(ncid, 'y', NF90_DOUBLE, &
            (/rho_dimid, theta_dimid, zeta_dimid/), y_varid)
        ierr = nf90_def_var(ncid, 'z', NF90_DOUBLE, &
            (/rho_dimid, theta_dimid, zeta_dimid/), z_varid)
        ierr = nf90_def_var(ncid, 'nfp', NF90_INT, nfp_varid)

        ierr = nf90_enddef(ncid)

        ! Write data
        ierr = nf90_put_var(ncid, rho_varid, rho)
        ierr = nf90_put_var(ncid, theta_varid, theta)
        ierr = nf90_put_var(ncid, zeta_varid, zeta)
        ierr = nf90_put_var(ncid, x_varid, x)
        ierr = nf90_put_var(ncid, y_varid, y)
        ierr = nf90_put_var(ncid, z_varid, z)
        ierr = nf90_put_var(ncid, nfp_varid, nfp)

        ierr = nf90_close(ncid)
    end subroutine write_chartmap_nc

end program vmec_to_chartmap
