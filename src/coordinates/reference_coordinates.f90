module reference_coordinates

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, make_vmec_coordinate_system, &
                                  make_geoflux_coordinate_system, make_chartmap_coordinate_system
    use netcdf, only: nf90_close, nf90_inq_dimid, nf90_inq_varid, nf90_noerr, nf90_nowrite, nf90_open

    implicit none

    class(coordinate_system_t), allocatable, public :: ref_coords

contains

    subroutine init_reference_coordinates(coord_input)
        character(*), intent(in) :: coord_input

        if (len_trim(coord_input) == 0) then
            print *, 'reference_coordinates.init_reference_coordinates: ', &
                'coord_input must be set (see params.apply_config_aliases)'
            error stop
        end if

        if (allocated(ref_coords)) deallocate (ref_coords)

        if (is_chartmap_file(coord_input)) then
            call make_chartmap_coordinate_system(ref_coords, trim(coord_input))
        else if (is_geqdsk_name(coord_input)) then
            call make_geoflux_coordinate_system(ref_coords)
        else
            call make_vmec_coordinate_system(ref_coords)
        end if
    end subroutine init_reference_coordinates

    logical function is_chartmap_file(filename)
        character(*), intent(in) :: filename

        integer :: ncid
        integer :: dimid
        integer :: varid
        integer :: ierr
        logical :: exists
        character(:), allocatable :: lower_name

        lower_name = to_lower(trim(filename))

        inquire(file=trim(filename), exist=exists)
        if (.not. exists) then
            is_chartmap_file = index(strip_directory(lower_name), 'chartmap') > 0
            return
        end if

        ierr = nf90_open(trim(filename), nf90_nowrite, ncid)
        if (ierr /= nf90_noerr) then
            is_chartmap_file = index(strip_directory(lower_name), 'chartmap') > 0
            return
        end if

        is_chartmap_file = .true.

        ierr = nf90_inq_dimid(ncid, 'rho', dimid)
        if (ierr /= nf90_noerr) is_chartmap_file = .false.
        ierr = nf90_inq_dimid(ncid, 'theta', dimid)
        if (ierr /= nf90_noerr) is_chartmap_file = .false.
        ierr = nf90_inq_dimid(ncid, 'zeta', dimid)
        if (ierr /= nf90_noerr) is_chartmap_file = .false.

        ierr = nf90_inq_varid(ncid, 'x', varid)
        if (ierr /= nf90_noerr) is_chartmap_file = .false.
        ierr = nf90_inq_varid(ncid, 'y', varid)
        if (ierr /= nf90_noerr) is_chartmap_file = .false.
        ierr = nf90_inq_varid(ncid, 'z', varid)
        if (ierr /= nf90_noerr) is_chartmap_file = .false.

        ierr = nf90_close(ncid)

        if (.not. is_chartmap_file) then
            is_chartmap_file = index(strip_directory(lower_name), 'chartmap') > 0
        end if
    end function is_chartmap_file

    logical function is_geqdsk_name(filename)
        character(*), intent(in) :: filename

        character(:), allocatable :: lower_name

        lower_name = to_lower(trim(filename))

        is_geqdsk_name = endswith(lower_name, '.geqdsk') .or. &
            endswith(lower_name, '.eqdsk')
        if (.not. is_geqdsk_name) then
            is_geqdsk_name = startswith(strip_directory(lower_name), 'geqdsk')
        end if
    end function is_geqdsk_name

    logical function startswith(text, start)
        character(*), intent(in) :: text
        character(*), intent(in) :: start
        integer :: len_text, len_start

        len_text = len_trim(text)
        len_start = len_trim(start)

        startswith = .false.
        if (len_text >= len_start) then
            startswith = (text(1:len_start) == start)
        end if
    end function startswith

    logical function endswith(text, ending)
        character(*), intent(in) :: text
        character(*), intent(in) :: ending
        integer :: len_text, len_end

        len_text = len_trim(text)
        len_end = len_trim(ending)

        endswith = .false.
        if (len_text >= len_end) then
            endswith = (text(len_text - len_end + 1:len_text) == ending)
        end if
    end function endswith

    function strip_directory(filename)
        character(*), intent(in) :: filename
        character(len(filename)) :: strip_directory
        integer :: i

        strip_directory = filename
        do i = len(filename), 1, -1
            if (filename(i:i) == '/') then
                strip_directory = filename(i + 1:len(filename))
                return
            end if
        end do
    end function strip_directory

    function to_lower(text) result(lower)
        character(*), intent(in) :: text
        character(len(text)) :: lower
        integer :: i

        lower = text
        do i = 1, len(text)
            select case (text(i:i))
            case ('A':'Z')
                lower(i:i) = achar(iachar(text(i:i)) + 32)
            end select
        end do
    end function to_lower

end module reference_coordinates
