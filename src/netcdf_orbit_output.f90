module netcdf_orbit_output
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use netcdf
    implicit none

    private
    public :: init_orbit_netcdf, write_orbit_step, flush_orbit, close_orbit_netcdf

    ! NetCDF file handle and variable IDs
    integer :: ncid
    integer :: dimid_particle, dimid_timestep
    integer :: varid_time, varid_s, varid_theta, varid_phi, varid_p_abs, varid_v_par
    logical :: netcdf_initialized = .false.
    real(dp) :: fill_value
    integer :: n_timesteps_max

    ! Per-thread orbit buffers
    real(dp), allocatable :: buffer_time(:,:)
    real(dp), allocatable :: buffer_s(:,:)
    real(dp), allocatable :: buffer_theta(:,:)
    real(dp), allocatable :: buffer_phi(:,:)
    real(dp), allocatable :: buffer_p_abs(:,:)
    real(dp), allocatable :: buffer_v_par(:,:)
    !$omp threadprivate(buffer_time, buffer_s, buffer_theta, buffer_phi, buffer_p_abs, buffer_v_par)

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
        status = nf90_put_att(ncid, varid_time, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att time _FillValue')

        status = nf90_def_var(ncid, 's', nf90_double, &
            [dimid_particle, dimid_timestep], varid_s)
        call check_nc(status, 'nf90_def_var s')
        status = nf90_put_att(ncid, varid_s, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att s _FillValue')

        status = nf90_def_var(ncid, 'theta', nf90_double, &
            [dimid_particle, dimid_timestep], varid_theta)
        call check_nc(status, 'nf90_def_var theta')
        status = nf90_put_att(ncid, varid_theta, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att theta _FillValue')

        status = nf90_def_var(ncid, 'phi', nf90_double, &
            [dimid_particle, dimid_timestep], varid_phi)
        call check_nc(status, 'nf90_def_var phi')
        status = nf90_put_att(ncid, varid_phi, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att phi _FillValue')

        status = nf90_def_var(ncid, 'p_abs', nf90_double, &
            [dimid_particle, dimid_timestep], varid_p_abs)
        call check_nc(status, 'nf90_def_var p_abs')
        status = nf90_put_att(ncid, varid_p_abs, '_FillValue', fill_value)
        call check_nc(status, 'nf90_put_att p_abs _FillValue')

        status = nf90_def_var(ncid, 'v_par', nf90_double, &
            [dimid_particle, dimid_timestep], varid_v_par)
        call check_nc(status, 'nf90_def_var v_par')
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

        ! Store max timesteps for buffer allocation
        n_timesteps_max = n_timesteps

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


    subroutine write_orbit_step(ipart, istep, t, xref, p_abs, v_par)
        integer, intent(in) :: ipart, istep
        real(dp), intent(in) :: t, xref(3), p_abs, v_par

        if (.not. netcdf_initialized) then
            error stop 'write_orbit_step: NetCDF not initialized'
        end if

        ! Allocate thread-local buffer on first use
        if (.not. allocated(buffer_time)) then
            allocate(buffer_time(1, n_timesteps_max))
            allocate(buffer_s(1, n_timesteps_max))
            allocate(buffer_theta(1, n_timesteps_max))
            allocate(buffer_phi(1, n_timesteps_max))
            allocate(buffer_p_abs(1, n_timesteps_max))
            allocate(buffer_v_par(1, n_timesteps_max))
            buffer_time(:,:) = fill_value
            buffer_s(:,:) = fill_value
            buffer_theta(:,:) = fill_value
            buffer_phi(:,:) = fill_value
            buffer_p_abs(:,:) = fill_value
            buffer_v_par(:,:) = fill_value
        end if

        ! Store data in thread-local buffer
        buffer_time(1, istep) = t
        buffer_s(1, istep) = xref(1)
        buffer_theta(1, istep) = xref(2)
        buffer_phi(1, istep) = xref(3)
        buffer_p_abs(1, istep) = p_abs
        buffer_v_par(1, istep) = v_par

    end subroutine write_orbit_step


    subroutine flush_orbit(ipart)
        integer, intent(in) :: ipart
        integer :: status
        integer :: start(2), count(2)

        if (.not. netcdf_initialized) then
            error stop 'flush_orbit: NetCDF not initialized'
        end if

        if (.not. allocated(buffer_time)) then
            return  ! No data to flush
        end if

        ! Write entire orbit at once
        start = [ipart, 1]
        count = [1, n_timesteps_max]

        !$omp critical (netcdf_write)

        status = nf90_put_var(ncid, varid_time, buffer_time, start=start, count=count)
        call check_nc(status, 'nf90_put_var time')

        status = nf90_put_var(ncid, varid_s, buffer_s, start=start, count=count)
        call check_nc(status, 'nf90_put_var s')

        status = nf90_put_var(ncid, varid_theta, buffer_theta, start=start, count=count)
        call check_nc(status, 'nf90_put_var theta')

        status = nf90_put_var(ncid, varid_phi, buffer_phi, start=start, count=count)
        call check_nc(status, 'nf90_put_var phi')

        status = nf90_put_var(ncid, varid_p_abs, buffer_p_abs, start=start, count=count)
        call check_nc(status, 'nf90_put_var p_abs')

        status = nf90_put_var(ncid, varid_v_par, buffer_v_par, start=start, count=count)
        call check_nc(status, 'nf90_put_var v_par')

        !$omp end critical (netcdf_write)

    end subroutine flush_orbit


    subroutine close_orbit_netcdf()
        integer :: status

        if (.not. netcdf_initialized) then
            return  ! Not an error, just nothing to close
        end if

        status = nf90_close(ncid)
        call check_nc(status, 'nf90_close')

        netcdf_initialized = .false.
    end subroutine close_orbit_netcdf

end module netcdf_orbit_output
