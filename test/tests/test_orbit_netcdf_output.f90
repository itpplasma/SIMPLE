program test_orbit_netcdf_output
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use netcdf
    implicit none

    character(len=*), parameter :: test_dir = '../../../examples/orbit_output_test'
    character(len=*), parameter :: nc_file = 'orbits.nc'
    integer :: status, ncid, dimid_particle, dimid_timestep
    integer :: n_particles, n_timesteps
    integer :: varid_time, varid_s, varid_theta, varid_phi, varid_p_abs, varid_v_par
    real(dp), allocatable :: time_data(:,:)
    real(dp), allocatable :: first_position(:, :)
    real(dp) :: start_state(5), angle_error(2), periods(2)
    integer :: field_ncid, field_varid, nfp
    integer :: start_unit, particle, component, position_ids(3)
    logical :: file_exists
    integer :: expected_particles = 32
    integer :: expected_timesteps = 10001

    print *, 'Testing NetCDF orbit output system...'

    ! Check if orbits.nc exists in test directory
    inquire(file=trim(test_dir)//'/'//trim(nc_file), exist=file_exists)
    if (.not. file_exists) then
        print *, 'ERROR: orbits.nc not found in ', trim(test_dir)
        error stop 1
    end if
    print *, 'PASS: orbits.nc file exists'

    ! Open NetCDF file
    status = nf90_open(trim(test_dir)//'/'//trim(nc_file), nf90_nowrite, ncid)
    if (status /= nf90_noerr) then
        print *, 'ERROR: Cannot open orbits.nc: ', trim(nf90_strerror(status))
        error stop 1
    end if
    print *, 'PASS: orbits.nc opened successfully'

    ! Check dimensions
    status = nf90_inq_dimid(ncid, 'particle', dimid_particle)
    if (status /= nf90_noerr) then
        print *, 'ERROR: particle dimension not found'
        status = nf90_close(ncid)
        error stop 1
    end if

    status = nf90_inquire_dimension(ncid, dimid_particle, len=n_particles)
    if (status /= nf90_noerr) then
        print *, 'ERROR: Cannot read particle dimension size'
        status = nf90_close(ncid)
        error stop 1
    end if

    if (n_particles /= expected_particles) then
        print *, 'ERROR: Expected ', expected_particles, ' particles, got ', n_particles
        status = nf90_close(ncid)
        error stop 1
    end if
    print *, 'PASS: particle dimension = ', n_particles

    status = nf90_inq_dimid(ncid, 'timestep', dimid_timestep)
    if (status /= nf90_noerr) then
        print *, 'ERROR: timestep dimension not found'
        status = nf90_close(ncid)
        error stop 1
    end if

    status = nf90_inquire_dimension(ncid, dimid_timestep, len=n_timesteps)
    if (status /= nf90_noerr) then
        print *, 'ERROR: Cannot read timestep dimension size'
        status = nf90_close(ncid)
        error stop 1
    end if

    if (n_timesteps /= expected_timesteps) then
        print *, 'ERROR: Expected ', expected_timesteps, ' timesteps, got ', n_timesteps
        status = nf90_close(ncid)
        error stop 1
    end if
    print *, 'PASS: timestep dimension = ', n_timesteps

    ! Check all required variables exist
    status = nf90_inq_varid(ncid, 'time', varid_time)
    if (status /= nf90_noerr) then
        print *, 'ERROR: time variable not found'
        status = nf90_close(ncid)
        error stop 1
    end if

    status = nf90_inq_varid(ncid, 's', varid_s)
    if (status /= nf90_noerr) then
        print *, 'ERROR: s variable not found'
        status = nf90_close(ncid)
        error stop 1
    end if

    status = nf90_inq_varid(ncid, 'theta', varid_theta)
    if (status /= nf90_noerr) then
        print *, 'ERROR: theta variable not found'
        status = nf90_close(ncid)
        error stop 1
    end if

    status = nf90_inq_varid(ncid, 'phi', varid_phi)
    if (status /= nf90_noerr) then
        print *, 'ERROR: phi variable not found'
        status = nf90_close(ncid)
        error stop 1
    end if

    status = nf90_inq_varid(ncid, 'p_abs', varid_p_abs)
    if (status /= nf90_noerr) then
        print *, 'ERROR: p_abs variable not found'
        status = nf90_close(ncid)
        error stop 1
    end if

    status = nf90_inq_varid(ncid, 'v_par', varid_v_par)
    if (status /= nf90_noerr) then
        print *, 'ERROR: v_par variable not found'
        status = nf90_close(ncid)
        error stop 1
    end if
    print *, 'PASS: All required variables exist'

    status = nf90_open(test_dir//'/wout.nc', nf90_nowrite, field_ncid)
    if (status /= nf90_noerr) error stop 'Cannot open orbit equilibrium'
    status = nf90_inq_varid(field_ncid, 'nfp', field_varid)
    if (status /= nf90_noerr) error stop 'Cannot find field-period count'
    status = nf90_get_var(field_ncid, field_varid, nfp)
    if (status /= nf90_noerr .or. nfp < 1) error stop 'Invalid field-period count'
    status = nf90_close(field_ncid)
    periods = 2.0_dp*acos(-1.0_dp)/[1.0_dp, real(nfp, dp)]

    ! The writer must recover the reference birth within each field period.
    allocate (first_position(n_particles, 3))
    position_ids = [varid_s, varid_theta, varid_phi]
    do component = 1, 3
        status = nf90_get_var(ncid, position_ids(component), &
                              first_position(:, component:component), &
                              start=[1, 1], count=[n_particles, 1])
        if (status /= nf90_noerr) error stop 'Cannot read first orbit position'
    end do
    open (newunit=start_unit, file=test_dir//'/start.dat', status='old')
    do particle = 1, n_particles
        read (start_unit, *) start_state
        angle_error = modulo(first_position(particle, 2:3) - &
                             start_state(2:3) + 0.5_dp*periods, periods) - &
                      0.5_dp*periods
        if (abs(first_position(particle, 1) - start_state(1)) > 1.0e-9_dp .or. &
            any(abs(angle_error) > 1.0e-9_dp)) then
            print *, 'Initial orbit sample differs from reference birth', particle
            print *, 'Expected:', start_state(1:3)
            print *, 'Actual:', first_position(particle, :)
            error stop 'Initial orbit sample differs from reference birth'
        end if
    end do
    close (start_unit)
    deallocate (first_position)

    ! Read a sample of time data to verify it's not all NaN
    allocate(time_data(n_particles, n_timesteps))
    status = nf90_get_var(ncid, varid_time, time_data)
    if (status /= nf90_noerr) then
        print *, 'ERROR: Cannot read time data'
        deallocate(time_data)
        status = nf90_close(ncid)
        error stop 1
    end if

    ! Check that some time data is valid (particle 3, timestep 1)
    ! Note: NetCDF stores as time(timestep, particle) but Fortran reads as time(particle, timestep)
    if (ieee_is_nan(time_data(3, 1))) then
        print *, 'ERROR: Time value for particle 3, timestep 1 is NaN'
        deallocate(time_data)
        status = nf90_close(ncid)
        error stop 1
    end if
    print *, 'PASS: Time data contains valid values'
    print *, '      Sample time value (particle 3, timestep 1): ', time_data(3, 1)

    deallocate(time_data)

    ! Close file
    status = nf90_close(ncid)
    if (status /= nf90_noerr) then
        print *, 'WARNING: Error closing file: ', trim(nf90_strerror(status))
    end if

    print *, ''
    print *, '========================================='
    print *, 'All NetCDF orbit output tests PASSED!'
    print *, '========================================='

end program test_orbit_netcdf_output
