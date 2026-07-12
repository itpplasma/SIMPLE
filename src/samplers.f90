module samplers
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util

    implicit none

    character(len=*), parameter :: START_FILE = 'start.dat'
    character(len=*), parameter :: START_FILE_ANTS = 'start_ants.dat'
    character(len=*), parameter :: START_FILE_BATCH = 'batch.dat'

    ! Interface ################################
    INTERFACE sample
        MODULE PROCEDURE sample_read
        MODULE PROCEDURE sample_surface_fieldline
        MODULE PROCEDURE sample_grid
        MODULE PROCEDURE sample_volume_single
        MODULE PROCEDURE sample_random_batch
        MODULE PROCEDURE sample_points_ants
    END INTERFACE sample

contains
    ! Functions #################################
    subroutine init_starting_surf
        use alpha_lifetime_sub, only: integrate_mfl_can
        use params, only: dphi, nper, npoiper, phibeg, thetabeg, volstart, &
                          xstart, sbeg, bmin, bmax, bmod00

        integer :: ierr = 0
        real(dp), dimension(npoiper*nper) :: bstart

        xstart = 0.d0
        bstart = 0.d0
        volstart = 0.d0

        ! For VMEC-backed runs the driver calls this while VMEC magfie is active,
        ! so xstart can be copied directly to reference-coordinate zstart. The
        ! volstart integral gives volume-weighted sampling on this one surface.
        call integrate_mfl_can( &
            npoiper*nper, dphi, sbeg(1), phibeg, thetabeg, &
            xstart, bstart, volstart, bmod00, ierr)

        if (ierr .ne. 0) then
            print *, 'starting field line has points outside the chamber'
            stop
        end if

        ! maximum value of B module:
        bmax = maxval(bstart)
        bmin = minval(bstart)

        print *, 'bmod00 = ', bmod00, 'bmin = ', bmin, 'bmax = ', bmax
    end subroutine init_starting_surf

    subroutine load_starting_points(zstart, filename)
        real(dp), dimension(:, :), intent(inout) :: zstart
        character(len=*), intent(in) :: filename
        integer :: ipart

        open (1, file=filename, recl=1024)
        do ipart = 1, size(zstart, 2)
            read (1, *) zstart(:, ipart)
        end do
        close (1)
    end subroutine load_starting_points

    subroutine save_starting_points(zstart)
        real(dp), dimension(:, :), intent(in) :: zstart
        integer :: ipart

        open (1, file=START_FILE, recl=1024)
        do ipart = 1, size(zstart, 2)
            write (1, *) zstart(:, ipart)
        end do
        close (1)
    end subroutine save_starting_points

    subroutine sample_read(zstart, filename)
        real(dp), dimension(:, :), intent(inout) :: zstart
        character(len=*), intent(in) :: filename

        call load_starting_points(zstart, filename)
    end subroutine

    ! Samplers ################################
    subroutine sample_volume_single(zstart, s_inner, s_outer)
        use params, only: isw_field_type, num_surf
        use field_can_mod, only: integ_to_ref

        real(dp), intent(in) :: s_inner
        real(dp), intent(in) :: s_outer
        real(dp), parameter :: s_min = 0.01d0
        real(dp) :: tmp_rand, s_lo, s_hi
        real(dp) :: r, vartheta, varphi
        real(dp), dimension(:, :), intent(inout) :: zstart
        integer :: ipart

        ! If user wants to do volume with 0 or 1 surfaces,
        !   we "add" the constraints, therefore having 2 surfaces.
        if (2 /= num_surf) then
            num_surf = 2
        end if

        ! Clamp lower bound to s_min to avoid axis singularity
        s_lo = max(s_inner, s_min)
        s_hi = max(s_outer, s_min)

        do ipart = 1, size(zstart, 2)
            call random_number(tmp_rand)
            r = tmp_rand*(s_hi - s_lo) + s_lo

            call random_number(tmp_rand)
            vartheta = twopi*tmp_rand
            call random_number(tmp_rand)
            varphi = twopi*tmp_rand
            ! we store starting points in reference coordinates:
            call integ_to_ref([r, vartheta, varphi], zstart(1:3, ipart))
            ! normalized velocity module z(4) = v / v_0:
            zstart(4, ipart) = 1.d0
            ! starting pitch z(5)=v_\parallel / v:
            call random_number(tmp_rand)
            zstart(5, ipart) = 2.d0*(tmp_rand - 0.5d0)
        end do

        call save_starting_points(zstart)

    end subroutine sample_volume_single

    subroutine sample_surface_fieldline(zstart)
        real(dp), dimension(:, :), intent(inout) :: zstart

        call sample_surface_fieldline_impl(zstart, .false.)
    end subroutine sample_surface_fieldline

    subroutine sample_surface_fieldline_from_integ(zstart)
        real(dp), dimension(:, :), intent(inout) :: zstart

        call sample_surface_fieldline_impl(zstart, .true.)
    end subroutine sample_surface_fieldline_from_integ

    subroutine sample_surface_fieldline_impl(zstart, xstart_is_integ_coords)
        use params, only: volstart, ibins, xstart, npoiper, nper
        use binsrc_sub, only: binsrc
        use field_can_mod, only: integ_to_ref

        real(dp), dimension(:, :), intent(inout) :: zstart
        logical, intent(in) :: xstart_is_integ_coords

        real(dp) :: xi
        integer :: ipart, i

        do ipart = 1, size(zstart, 2)
            call random_number(xi)
            call binsrc(volstart, 1, npoiper*nper, xi, i)
            ibins = i
            if (xstart_is_integ_coords) then
                call integ_to_ref(xstart(:, i), zstart(1:3, ipart))
            else
                zstart(1:3, ipart) = xstart(:, i)
            end if
            zstart(4, ipart) = 1.d0  ! normalized velocity module z(4) = v / v_0
            call random_number(xi)
            zstart(5, ipart) = 2.d0*(xi - 0.5d0)  ! starting pitch z(5)=v_\parallel / v
        end do

        call save_starting_points(zstart)

    end subroutine sample_surface_fieldline_impl

    subroutine sample_grid(zstart, grid_density, xstart_is_integ_coords)
        use params, only: ntestpart, zstart_dim1, zend, times_lost, &
                          orbit_exit_code, boundary_event_radial_residual, &
                          boundary_event_time_width, trap_par, perp_inv, iclass, sbeg
        use util, only: pi
        use field_can_mod, only: integ_to_ref

        real(dp), dimension(:, :), allocatable, intent(inout) :: zstart
        real(dp), intent(in) :: grid_density
        logical, intent(in), optional :: xstart_is_integ_coords
        real(dp) :: xi, xsize_real
        real(dp) :: xinteg(3)
        integer :: ngrid, ipart, jpart, lidx
        logical :: convert_surface_starts

        convert_surface_starts = .false.
        if (present(xstart_is_integ_coords)) then
            convert_surface_starts = xstart_is_integ_coords
        end if

        xsize_real = (2*pi)*grid_density !angle density
        ngrid = int((1/grid_density) - 1)
        ntestpart = ngrid**2 !number of total angle points

        ! Resize particle coord. arrays and result memory.
        if (allocated(zstart)) deallocate (zstart)
        if (allocated(zend)) deallocate (zend)
        allocate (zstart(zstart_dim1, ntestpart), zend(zstart_dim1, ntestpart))
        if (allocated(times_lost)) deallocate (times_lost)
        if (allocated(orbit_exit_code)) deallocate (orbit_exit_code)
        if (allocated(boundary_event_radial_residual)) &
            deallocate (boundary_event_radial_residual)
        if (allocated(boundary_event_time_width)) deallocate (boundary_event_time_width)
        if (allocated(trap_par)) deallocate (trap_par)
        if (allocated(perp_inv)) deallocate (perp_inv)
        if (allocated(iclass)) deallocate (iclass)
        allocate(times_lost(ntestpart), orbit_exit_code(ntestpart), &
                 boundary_event_radial_residual(ntestpart), &
                 boundary_event_time_width(ntestpart), &
                 trap_par(ntestpart), perp_inv(ntestpart), iclass(3,ntestpart))
        orbit_exit_code = 0
        boundary_event_radial_residual = -1d0
        boundary_event_time_width = -1d0

        do ipart = 1, ngrid
            do jpart = 1, ngrid
                lidx = (jpart - 1)*ngrid + ipart
                xinteg = [sbeg(1), xsize_real*ipart, xsize_real*jpart]
                if (convert_surface_starts) then
                    call integ_to_ref(xinteg, zstart(1:3, lidx))
                else
                    zstart(1:3, lidx) = xinteg
                end if
                zstart(4, lidx) = 1.d0  ! normalized velocity module z(4) = v / v_0
                call random_number(xi)
                zstart(5, lidx) = 2.d0*(xi - 0.5d0)  ! starting pitch z(5)=v_\parallel / v
            end do
        end do

        call save_starting_points(zstart)

    end subroutine sample_grid

    subroutine sample_random_batch(zstart, reuse_existing)
        ! Get random batch from preexisting zstart, allows reuse.
        use params, only: batch_size, ntestpart, zstart_dim1, idx

        integer :: ran_begin, ran_end, ipart
        real :: temp_ran
        real(dp), dimension(:, :), intent(inout) :: zstart
        real(dp), dimension(zstart_dim1, batch_size) :: zstart_batch
        logical, intent(in) :: reuse_existing

        if (reuse_existing .eqv. .True.) then
            call load_starting_points(zstart_batch, START_FILE_BATCH)
        else
            call load_starting_points(zstart_batch, START_FILE)
            call random_number(temp_ran)
            ran_begin = INT(temp_ran)
            ran_end = ran_begin + batch_size
            if ((ran_end) .gt. (ntestpart)) then
                ran_begin = ran_begin - (ran_end - ntestpart)
            end if
            do ipart = 0, batch_size
                zstart(:, ipart) = zstart_batch(:, (ipart + ran_begin))
            end do
        end if

        do ipart = idx(0), idx(ntestpart)
            read (1, *) zstart(:, ipart)
        end do

    end subroutine sample_random_batch

    subroutine sample_points_ants(use_special_ants_file)
        use parse_ants, only: process_line
        use get_can_sub, only: vmec_to_can
        use params, only: ntestpart, zstart ! ANTS sampler uses global zstart

        logical, intent(in) :: use_special_ants_file

        integer, parameter :: maxlen = 4096
        character(len=maxlen) :: line
        real(8) :: v_par, v_perp, u, v, s
        real(8) :: th, ph
        integer :: ipart

        do ipart = 1, ntestpart
            if (use_special_ants_file) then
                open (1, file=START_FILE_ANTS, recl=1024)
                read (1, '(A)') line
                close (1)
            else
                open (1, file=START_FILE, recl=1024)
                read (1, '(A)') line
                close (1)
            end if

            call process_line(line, v_par, v_perp, u, v, s)
            ! In the test case, u runs from 0 to 1 and v from 0 to 4
            th = 2d0*pi*u
            ph = 2d0*pi*v/4d0
            zstart(1, ipart) = s
            zstart(2, ipart) = th
            zstart(3, ipart) = ph
            zstart(4, ipart) = 1.d0
            zstart(5, ipart) = v_par/sqrt(v_par**2 + v_perp**2)
        end do
    end subroutine sample_points_ants

end module samplers
