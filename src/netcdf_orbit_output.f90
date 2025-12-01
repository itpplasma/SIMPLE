module netcdf_orbit_output
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use netcdf
    implicit none

    private
    public :: init_orbit_netcdf, close_orbit_netcdf
    public :: write_orbit_to_netcdf
    public :: netcdf_compression_level

    ! Compression level: 0=none, 1-9=deflate (4 recommended)
    integer :: netcdf_compression_level = 4

    ! NetCDF file handle and variable IDs
    integer :: ncid
    integer :: dimid_particle, dimid_timestep
    integer :: varid_time, varid_s, varid_theta, varid_phi, varid_p_abs, varid_v_par
    logical :: netcdf_initialized = .false.
    real(dp) :: fill_value

contains

    subroutine check_nc(status, location)
        integer, intent(in) :: status
        character(len=*), intent(in) :: location

        if (status /= nf90_noerr) then
            print *, 'NetCDF error at ', trim(location), ': ', trim(nf90_strerror(status))
            error stop 'NetCDF operation failed'
        end if
    end subroutine check_nc


    subroutine init_orbit_netcdf(n_particles, n_timesteps)
        integer, intent(in) :: n_particles, n_timesteps
        integer :: status, coord_varid_particle, coord_varid_timestep

        if (netcdf_initialized) then
            error stop 'init_orbit_netcdf: already initialized'
        end if

        fill_value = ieee_value(0.0d0, ieee_quiet_nan)

        ! Create NetCDF file
        status = nf90_create('orbits.nc', nf90_netcdf4, ncid)
        call check_nc(status, 'nf90_create')

        ! Define dimensions
        status = nf90_def_dim(ncid, 'particle', n_particles, dimid_particle)
        call check_nc(status, 'nf90_def_dim particle')

        status = nf90_def_dim(ncid, 'timestep', n_timesteps, dimid_timestep)
        call check_nc(status, 'nf90_def_dim timestep')

        ! Define coordinate variables
        status = nf90_def_var(ncid, 'particle', nf90_int, [dimid_particle], coord_varid_particle)
        call check_nc(status, 'nf90_def_var particle')

        status = nf90_def_var(ncid, 'timestep', nf90_int, [dimid_timestep], coord_varid_timestep)
        call check_nc(status, 'nf90_def_var timestep')

        ! Define data variables (all 2D: particle x timestep)
        status = nf90_def_var(ncid, 'time', nf90_double, &
            [dimid_particle, dimid_timestep], varid_time)
        call check_nc(status, 'nf90_def_var time')
        if (netcdf_compression_level > 0) then
            status = nf90_def_var_deflate(ncid, varid_time, 1, 1, netcdf_compression_level)
            call check_nc(status, 'nf90_def_var_deflate time')
        end if
        status = nf90_put_att(ncid, varid_time, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att time _FillValue')

        status = nf90_def_var(ncid, 's', nf90_double, &
            [dimid_particle, dimid_timestep], varid_s)
        call check_nc(status, 'nf90_def_var s')
        if (netcdf_compression_level > 0) then
            status = nf90_def_var_deflate(ncid, varid_s, 1, 1, netcdf_compression_level)
            call check_nc(status, 'nf90_def_var_deflate s')
        end if
        status = nf90_put_att(ncid, varid_s, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att s _FillValue')

        status = nf90_def_var(ncid, 'theta', nf90_double, &
            [dimid_particle, dimid_timestep], varid_theta)
        call check_nc(status, 'nf90_def_var theta')
        if (netcdf_compression_level > 0) then
            status = nf90_def_var_deflate(ncid, varid_theta, 1, 1, netcdf_compression_level)
            call check_nc(status, 'nf90_def_var_deflate theta')
        end if
        status = nf90_put_att(ncid, varid_theta, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att theta _FillValue')

        status = nf90_def_var(ncid, 'phi', nf90_double, &
            [dimid_particle, dimid_timestep], varid_phi)
        call check_nc(status, 'nf90_def_var phi')
        if (netcdf_compression_level > 0) then
            status = nf90_def_var_deflate(ncid, varid_phi, 1, 1, netcdf_compression_level)
            call check_nc(status, 'nf90_def_var_deflate phi')
        end if
        status = nf90_put_att(ncid, varid_phi, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att phi _FillValue')

        status = nf90_def_var(ncid, 'p_abs', nf90_double, &
            [dimid_particle, dimid_timestep], varid_p_abs)
        call check_nc(status, 'nf90_def_var p_abs')
        if (netcdf_compression_level > 0) then
            status = nf90_def_var_deflate(ncid, varid_p_abs, 1, 1, netcdf_compression_level)
            call check_nc(status, 'nf90_def_var_deflate p_abs')
        end if
        status = nf90_put_att(ncid, varid_p_abs, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att p_abs _FillValue')

        status = nf90_def_var(ncid, 'v_par', nf90_double, &
            [dimid_particle, dimid_timestep], varid_v_par)
        call check_nc(status, 'nf90_def_var v_par')
        if (netcdf_compression_level > 0) then
            status = nf90_def_var_deflate(ncid, varid_v_par, 1, 1, netcdf_compression_level)
            call check_nc(status, 'nf90_def_var_deflate v_par')
        end if
        status = nf90_put_att(ncid, varid_v_par, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att v_par _FillValue')

        ! Add global attributes
        status = nf90_put_att(ncid, nf90_global, 'description', &
            'SIMPLE orbit tracer output (macrostep mode)')
        call check_nc(status, 'nf90_put_att description')

        status = nf90_put_att(ncid, nf90_global, 'source', &
            'Direct Fortran NetCDF output')
        call check_nc(status, 'nf90_put_att source')

        status = nf90_put_att(ncid, nf90_global, 'n_particles', n_particles)
        call check_nc(status, 'nf90_put_att n_particles')

        status = nf90_put_att(ncid, nf90_global, 'n_timesteps', n_timesteps)
        call check_nc(status, 'nf90_put_att n_timesteps')

        ! End define mode
        status = nf90_enddef(ncid)
        call check_nc(status, 'nf90_enddef')

        ! Write coordinate variables
        call write_coordinates(n_particles, n_timesteps)

        netcdf_initialized = .true.
    end subroutine init_orbit_netcdf


    subroutine write_coordinates(n_particles, n_timesteps)
        integer, intent(in) :: n_particles, n_timesteps
        integer :: status, i
        integer :: coord_varid_particle, coord_varid_timestep
        integer, allocatable :: particle_ids(:), timestep_ids(:)

        allocate(particle_ids(n_particles))
        allocate(timestep_ids(n_timesteps))

        do i = 1, n_particles
            particle_ids(i) = i
        end do

        do i = 1, n_timesteps
            timestep_ids(i) = i - 1  ! 0-based indexing for timesteps
        end do

        status = nf90_inq_varid(ncid, 'particle', coord_varid_particle)
        call check_nc(status, 'nf90_inq_varid particle')
        status = nf90_put_var(ncid, coord_varid_particle, particle_ids)
        call check_nc(status, 'nf90_put_var particle')

        status = nf90_inq_varid(ncid, 'timestep', coord_varid_timestep)
        call check_nc(status, 'nf90_inq_varid timestep')
        status = nf90_put_var(ncid, coord_varid_timestep, timestep_ids)
        call check_nc(status, 'nf90_put_var timestep')

        deallocate(particle_ids, timestep_ids)
    end subroutine write_coordinates




    subroutine write_orbit_to_netcdf(ipart, orbit_traj, orbit_times)
        use field_can_mod, only : can_to_ref
        use velo_mod, only : isw_field_type
        use magfie_sub, only : COILS
        use simple_coordinates, only : transform_cart_to_vmec
        integer, intent(in) :: ipart
        real(kind(1.0d0)), intent(in) :: orbit_traj(:,:)  ! (5, ntimstep)
        real(kind(1.0d0)), intent(in) :: orbit_times(:)   ! (ntimstep)

        integer :: it, status, n_times
        real(kind(1.0d0)) :: xref(3)
        real(kind(1.0d0)), allocatable :: s_array(:), theta_array(:), phi_array(:)

        if (.not. netcdf_initialized) then
            error stop 'write_orbit_to_netcdf: NetCDF not initialized'
        end if

        n_times = size(orbit_times)
        allocate(s_array(n_times), theta_array(n_times), phi_array(n_times))

        ! Convert canonical to reference coordinates
        do it = 1, n_times
            if (isw_field_type == COILS) then
                call transform_cart_to_vmec(orbit_traj(1:3, it), xref)
            else
                call can_to_ref(orbit_traj(1:3, it), xref)
            end if
            s_array(it) = xref(1)
            theta_array(it) = xref(2)
            phi_array(it) = xref(3)
        enddo

        ! Write entire orbit to NetCDF in one shot
        status = nf90_put_var(ncid, varid_time, orbit_times, start=[ipart, 1], count=[1, n_times])
        call check_nc(status, 'nf90_put_var time')

        status = nf90_put_var(ncid, varid_s, s_array, start=[ipart, 1], count=[1, n_times])
        call check_nc(status, 'nf90_put_var s')

        status = nf90_put_var(ncid, varid_theta, theta_array, start=[ipart, 1], count=[1, n_times])
        call check_nc(status, 'nf90_put_var theta')

        status = nf90_put_var(ncid, varid_phi, phi_array, start=[ipart, 1], count=[1, n_times])
        call check_nc(status, 'nf90_put_var phi')

        status = nf90_put_var(ncid, varid_p_abs, orbit_traj(4, :), start=[ipart, 1], count=[1, n_times])
        call check_nc(status, 'nf90_put_var p_abs')

        status = nf90_put_var(ncid, varid_v_par, orbit_traj(5, :), start=[ipart, 1], count=[1, n_times])
        call check_nc(status, 'nf90_put_var v_par')
    end subroutine write_orbit_to_netcdf


    subroutine close_orbit_netcdf()
        integer :: status

        if (.not. netcdf_initialized) then
            return  ! Not an error, just nothing to close
        end if

        ! Data already written incrementally, just close the file
        status = nf90_close(ncid)
        call check_nc(status, 'nf90_close')

        netcdf_initialized = .false.
    end subroutine close_orbit_netcdf

end module netcdf_orbit_output
