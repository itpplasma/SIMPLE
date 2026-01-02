module netcdf_results_output
    !> NetCDF output for particle tracing results.
    !>
    !> Provides structured binary output (results.nc) containing:
    !>   - Per-particle loss data (times_lost, trap_par, perp_inv, classification)
    !>   - Phase space positions (zstart, zend) in reference coordinates
    !>   - Cartesian positions (xstart_cart, xend_cart) for visualization/heat loads
    !>   - Simulation config as global attributes
    !>
    !> Usage:
    !>   call write_results_netcdf(filename)
    !>
    !> The Cartesian positions decouple post-processing from coordinate
    !> transformation machinery - Python can directly read (x,y,z) positions.

    use, intrinsic :: iso_fortran_env, only: dp => real64, int8
    use netcdf
    implicit none

    private
    public :: write_results_netcdf

contains

    subroutine check_nc(status, location)
        integer, intent(in) :: status
        character(len=*), intent(in) :: location

        if (status /= nf90_noerr) then
            print *, 'NetCDF error at ', trim(location), ': ', &
                trim(nf90_strerror(status))
            error stop 'NetCDF operation failed'
        end if
    end subroutine check_nc

    subroutine write_results_netcdf(filename)
        !> Write all particle results to NetCDF file.
        !>
        !> Computes Cartesian positions from reference coordinates and writes
        !> complete results including config attributes.

        use params, only: ntestpart, trace_time, zstart, zend, times_lost, &
                          trap_par, perp_inv, iclass, class_lost, &
                          facE_al, v0, integmode, relerr, npoiper, npoiper2, &
                          ntimstep, startmode, num_surf, sbeg, &
                          isw_field_type, swcoll, deterministic, ran_seed, &
                          netcdffile, field_input, coord_input
        use reference_coordinates, only: ref_coords

        character(len=*), intent(in) :: filename

        integer :: ncid, status
        integer :: dim_particle, dim_xyz, dim_phase, dim_iclass
        integer :: var_times_lost, var_trap_par, var_perp_inv
        integer :: var_zstart, var_zend, var_xstart_cart, var_xend_cart
        integer :: var_iclass, var_class_lost
        integer :: i
        real(dp) :: xstart_cart(3, ntestpart), xend_cart(3, ntestpart)
        integer(int8) :: class_lost_i8(ntestpart)
        character(len=32) :: coord_type

        ! Compute Cartesian positions
        call compute_cartesian_positions(xstart_cart, xend_cart)

        ! Convert logical to int8
        do i = 1, ntestpart
            if (class_lost(i)) then
                class_lost_i8(i) = 1_int8
            else
                class_lost_i8(i) = 0_int8
            end if
        end do

        ! Determine coordinate type
        if (len_trim(coord_input) > 0) then
            if (index(coord_input, 'chartmap') > 0 .or. &
                index(coord_input, '.nc') > 0) then
                coord_type = 'chartmap'
            else
                coord_type = 'vmec'
            end if
        else
            coord_type = 'vmec'
        end if

        ! Create file
        status = nf90_create(filename, nf90_netcdf4, ncid)
        call check_nc(status, 'nf90_create')

        ! Define dimensions
        status = nf90_def_dim(ncid, 'particle', ntestpart, dim_particle)
        call check_nc(status, 'def_dim particle')

        status = nf90_def_dim(ncid, 'xyz', 3, dim_xyz)
        call check_nc(status, 'def_dim xyz')

        status = nf90_def_dim(ncid, 'phase', 5, dim_phase)
        call check_nc(status, 'def_dim phase')

        status = nf90_def_dim(ncid, 'iclass_dim', 3, dim_iclass)
        call check_nc(status, 'def_dim iclass_dim')

        ! Define variables with full attribute error checking
        status = nf90_def_var(ncid, 'times_lost', nf90_double, &
            [dim_particle], var_times_lost)
        call check_nc(status, 'def_var times_lost')
        status = nf90_put_att(ncid, var_times_lost, 'units', 's')
        call check_nc(status, 'put_att times_lost units')
        status = nf90_put_att(ncid, var_times_lost, 'description', &
            'Loss time (-1 if confined)')
        call check_nc(status, 'put_att times_lost description')

        status = nf90_def_var(ncid, 'trap_par', nf90_double, &
            [dim_particle], var_trap_par)
        call check_nc(status, 'def_var trap_par')
        status = nf90_put_att(ncid, var_trap_par, 'units', 'dimensionless')
        call check_nc(status, 'put_att trap_par units')
        status = nf90_put_att(ncid, var_trap_par, 'description', &
            'Trapping parameter eta = mu*B_max/E')
        call check_nc(status, 'put_att trap_par description')

        status = nf90_def_var(ncid, 'perp_inv', nf90_double, &
            [dim_particle], var_perp_inv)
        call check_nc(status, 'def_var perp_inv')
        status = nf90_put_att(ncid, var_perp_inv, 'units', 'cm^2*G')
        call check_nc(status, 'put_att perp_inv units')
        status = nf90_put_att(ncid, var_perp_inv, 'description', &
            'Perpendicular adiabatic invariant J_perp = mu*B')
        call check_nc(status, 'put_att perp_inv description')

        status = nf90_def_var(ncid, 'zstart', nf90_double, &
            [dim_phase, dim_particle], var_zstart)
        call check_nc(status, 'def_var zstart')
        status = nf90_put_att(ncid, var_zstart, 'description', &
            'Initial phase space in reference coordinates')
        call check_nc(status, 'put_att zstart description')
        status = nf90_put_att(ncid, var_zstart, 'components', &
            's (flux label), theta (rad), zeta (rad), p (v/v0), v_par/v')
        call check_nc(status, 'put_att zstart components')
        status = nf90_put_att(ncid, var_zstart, 'coordinate_note', &
            'See coordinate_type global attribute for VMEC vs chartmap')
        call check_nc(status, 'put_att zstart coordinate_note')

        status = nf90_def_var(ncid, 'zend', nf90_double, &
            [dim_phase, dim_particle], var_zend)
        call check_nc(status, 'def_var zend')
        status = nf90_put_att(ncid, var_zend, 'description', &
            'Final phase space in reference coordinates (0 if not traced)')
        call check_nc(status, 'put_att zend description')
        status = nf90_put_att(ncid, var_zend, 'components', &
            's (flux label), theta (rad), zeta (rad), p (v/v0), v_par/v')
        call check_nc(status, 'put_att zend components')

        status = nf90_def_var(ncid, 'xstart_cart', nf90_double, &
            [dim_xyz, dim_particle], var_xstart_cart)
        call check_nc(status, 'def_var xstart_cart')
        status = nf90_put_att(ncid, var_xstart_cart, 'units', 'cm')
        call check_nc(status, 'put_att xstart_cart units')
        status = nf90_put_att(ncid, var_xstart_cart, 'description', &
            'Initial Cartesian position (x, y, z)')
        call check_nc(status, 'put_att xstart_cart description')

        status = nf90_def_var(ncid, 'xend_cart', nf90_double, &
            [dim_xyz, dim_particle], var_xend_cart)
        call check_nc(status, 'def_var xend_cart')
        status = nf90_put_att(ncid, var_xend_cart, 'units', 'cm')
        call check_nc(status, 'put_att xend_cart units')
        status = nf90_put_att(ncid, var_xend_cart, 'description', &
            'Final Cartesian position (equals xstart if not traced)')
        call check_nc(status, 'put_att xend_cart description')

        status = nf90_def_var(ncid, 'iclass', nf90_int, &
            [dim_iclass, dim_particle], var_iclass)
        call check_nc(status, 'def_var iclass')
        status = nf90_put_att(ncid, var_iclass, 'description', &
            'Classification: [J_par, topological, fractal] (0=unclassified)')
        call check_nc(status, 'put_att iclass description')

        status = nf90_def_var(ncid, 'class_lost', nf90_byte, &
            [dim_particle], var_class_lost)
        call check_nc(status, 'def_var class_lost')
        status = nf90_put_att(ncid, var_class_lost, 'description', &
            'Lost flag (1=lost, 0=confined)')
        call check_nc(status, 'put_att class_lost description')

        ! Global attributes - simulation config
        status = nf90_put_att(ncid, nf90_global, 'title', &
            'SIMPLE particle tracing results')
        call check_nc(status, 'put_att title')
        status = nf90_put_att(ncid, nf90_global, 'ntestpart', ntestpart)
        call check_nc(status, 'put_att ntestpart')
        status = nf90_put_att(ncid, nf90_global, 'trace_time', trace_time)
        call check_nc(status, 'put_att trace_time')
        status = nf90_put_att(ncid, nf90_global, 'facE_al', facE_al)
        call check_nc(status, 'put_att facE_al')
        status = nf90_put_att(ncid, nf90_global, 'v0', v0)
        call check_nc(status, 'put_att v0')
        status = nf90_put_att(ncid, nf90_global, 'integmode', integmode)
        call check_nc(status, 'put_att integmode')
        status = nf90_put_att(ncid, nf90_global, 'relerr', relerr)
        call check_nc(status, 'put_att relerr')
        status = nf90_put_att(ncid, nf90_global, 'npoiper', npoiper)
        call check_nc(status, 'put_att npoiper')
        status = nf90_put_att(ncid, nf90_global, 'npoiper2', npoiper2)
        call check_nc(status, 'put_att npoiper2')
        status = nf90_put_att(ncid, nf90_global, 'ntimstep', ntimstep)
        call check_nc(status, 'put_att ntimstep')
        status = nf90_put_att(ncid, nf90_global, 'startmode', startmode)
        call check_nc(status, 'put_att startmode')
        status = nf90_put_att(ncid, nf90_global, 'num_surf', num_surf)
        call check_nc(status, 'put_att num_surf')
        status = nf90_put_att(ncid, nf90_global, 'sbeg', sbeg)
        call check_nc(status, 'put_att sbeg')
        status = nf90_put_att(ncid, nf90_global, 'isw_field_type', isw_field_type)
        call check_nc(status, 'put_att isw_field_type')
        status = nf90_put_att(ncid, nf90_global, 'swcoll', merge(1, 0, swcoll))
        call check_nc(status, 'put_att swcoll')
        status = nf90_put_att(ncid, nf90_global, 'deterministic', &
            merge(1, 0, deterministic))
        call check_nc(status, 'put_att deterministic')
        status = nf90_put_att(ncid, nf90_global, 'ran_seed', ran_seed)
        call check_nc(status, 'put_att ran_seed')
        status = nf90_put_att(ncid, nf90_global, 'netcdffile', trim(netcdffile))
        call check_nc(status, 'put_att netcdffile')
        status = nf90_put_att(ncid, nf90_global, 'field_input', trim(field_input))
        call check_nc(status, 'put_att field_input')
        status = nf90_put_att(ncid, nf90_global, 'coord_input', trim(coord_input))
        call check_nc(status, 'put_att coord_input')
        status = nf90_put_att(ncid, nf90_global, 'coordinate_type', trim(coord_type))
        call check_nc(status, 'put_att coordinate_type')

        ! End define mode
        status = nf90_enddef(ncid)
        call check_nc(status, 'enddef')

        ! Write data
        status = nf90_put_var(ncid, var_times_lost, times_lost)
        call check_nc(status, 'put_var times_lost')

        status = nf90_put_var(ncid, var_trap_par, trap_par)
        call check_nc(status, 'put_var trap_par')

        status = nf90_put_var(ncid, var_perp_inv, perp_inv)
        call check_nc(status, 'put_var perp_inv')

        status = nf90_put_var(ncid, var_zstart, zstart)
        call check_nc(status, 'put_var zstart')

        status = nf90_put_var(ncid, var_zend, zend)
        call check_nc(status, 'put_var zend')

        status = nf90_put_var(ncid, var_xstart_cart, xstart_cart)
        call check_nc(status, 'put_var xstart_cart')

        status = nf90_put_var(ncid, var_xend_cart, xend_cart)
        call check_nc(status, 'put_var xend_cart')

        status = nf90_put_var(ncid, var_iclass, iclass)
        call check_nc(status, 'put_var iclass')

        status = nf90_put_var(ncid, var_class_lost, class_lost_i8)
        call check_nc(status, 'put_var class_lost')

        ! Close file
        status = nf90_close(ncid)
        call check_nc(status, 'close')

        print *, 'INFO: Results written to ', trim(filename)

    end subroutine write_results_netcdf

    subroutine compute_cartesian_positions(xstart_cart, xend_cart)
        !> Convert reference coordinates to Cartesian for all particles.
        !>
        !> For particles that were not traced (zend = 0), xend_cart is set
        !> to xstart_cart since they effectively remained at the start.
        !>
        !> TODO: Investigate batch coordinate transformation for performance.
        !> Current implementation uses O(n) individual calls.

        use params, only: ntestpart, zstart, zend
        use reference_coordinates, only: ref_coords

        real(dp), intent(out) :: xstart_cart(3, ntestpart)
        real(dp), intent(out) :: xend_cart(3, ntestpart)
        integer :: i
        logical :: zend_is_zero

        if (.not. allocated(ref_coords)) then
            print *, 'ERROR: ref_coords not allocated - cannot compute Cartesian'
            print *, 'This indicates init_reference_coordinates was not called'
            error stop 'NetCDF results output requires initialized coordinate system'
        end if

        do i = 1, ntestpart
            call ref_coords%evaluate_cart(zstart(1:3, i), xstart_cart(:, i))

            ! Check if particle was traced: zend is set to exactly 0 for untraced
            ! particles (see simple_main.f90:493 and classification.f90:94)
            zend_is_zero = all(zend(1:3, i) == 0.0_dp)
            if (zend_is_zero) then
                ! Particle not traced - use start position
                xend_cart(:, i) = xstart_cart(:, i)
            else
                call ref_coords%evaluate_cart(zend(1:3, i), xend_cart(:, i))
            end if
        end do

    end subroutine compute_cartesian_positions

end module netcdf_results_output
