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
    integer :: n_particles_total

    ! Shared buffer for all particles (write in one shot)
    real(dp), allocatable :: shared_buffer_time(:,:)
    real(dp), allocatable :: shared_buffer_s(:,:)
    real(dp), allocatable :: shared_buffer_theta(:,:)
    real(dp), allocatable :: shared_buffer_phi(:,:)
    real(dp), allocatable :: shared_buffer_p_abs(:,:)
    real(dp), allocatable :: shared_buffer_v_par(:,:)

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

        ! Store dimensions and allocate shared buffer
        n_timesteps_max = n_timesteps
        n_particles_total = n_particles

        allocate(shared_buffer_time(n_particles, n_timesteps))
        allocate(shared_buffer_s(n_particles, n_timesteps))
        allocate(shared_buffer_theta(n_particles, n_timesteps))
        allocate(shared_buffer_phi(n_particles, n_timesteps))
        allocate(shared_buffer_p_abs(n_particles, n_timesteps))
        allocate(shared_buffer_v_par(n_particles, n_timesteps))

        shared_buffer_time(:,:) = fill_value
        shared_buffer_s(:,:) = fill_value
        shared_buffer_theta(:,:) = fill_value
        shared_buffer_phi(:,:) = fill_value
        shared_buffer_p_abs(:,:) = fill_value
        shared_buffer_v_par(:,:) = fill_value

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

        ! Write directly to shared buffer (no critical section needed for different particles)
        shared_buffer_time(ipart, istep) = t
        shared_buffer_s(ipart, istep) = xref(1)
        shared_buffer_theta(ipart, istep) = xref(2)
        shared_buffer_phi(ipart, istep) = xref(3)
        shared_buffer_p_abs(ipart, istep) = p_abs
        shared_buffer_v_par(ipart, istep) = v_par

    end subroutine write_orbit_step


    subroutine flush_orbit(ipart)
        integer, intent(in) :: ipart
        ! No-op: data is already in shared buffer, will be written in close_orbit_netcdf
    end subroutine flush_orbit


    subroutine close_orbit_netcdf()
        use timing, only: get_wtime
        integer :: status
        real(dp) :: t_start, t_end

        if (.not. netcdf_initialized) then
            return  ! Not an error, just nothing to close
        end if

        ! Write entire buffer in one massive operation
        t_start = get_wtime()
        print *, 'Writing all orbit data to NetCDF...'

        status = nf90_put_var(ncid, varid_time, shared_buffer_time)
        call check_nc(status, 'nf90_put_var time (full)')

        status = nf90_put_var(ncid, varid_s, shared_buffer_s)
        call check_nc(status, 'nf90_put_var s (full)')

        status = nf90_put_var(ncid, varid_theta, shared_buffer_theta)
        call check_nc(status, 'nf90_put_var theta (full)')

        status = nf90_put_var(ncid, varid_phi, shared_buffer_phi)
        call check_nc(status, 'nf90_put_var phi (full)')

        status = nf90_put_var(ncid, varid_p_abs, shared_buffer_p_abs)
        call check_nc(status, 'nf90_put_var p_abs (full)')

        status = nf90_put_var(ncid, varid_v_par, shared_buffer_v_par)
        call check_nc(status, 'nf90_put_var v_par (full)')

        t_end = get_wtime()
        print '(A,F8.3,A)', 'NetCDF write completed in ', t_end - t_start, ' s'

        status = nf90_close(ncid)
        call check_nc(status, 'nf90_close')

        netcdf_initialized = .false.
    end subroutine close_orbit_netcdf

end module netcdf_orbit_output
