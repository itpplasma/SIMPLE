module field
    !> Field module aggregating all field types and factory functions.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: detect_refcoords_file_type, refcoords_file_chartmap, &
                                  refcoords_file_vmec_wout, refcoords_file_unknown
    use field_base, only: magnetic_field_t
    use field_vmec, only: vmec_field_t, create_vmec_field
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
#ifdef GVEC_AVAILABLE
    use field_gvec, only: gvec_field_t, create_gvec_field
#endif

    implicit none

contains

    subroutine field_from_file(filename, field)
        !> Create appropriate field type from file.
        !> For coils files, creates a splined_field_t wrapping coils_field_t.
        use reference_coordinates, only: ref_coords

        character(*), intent(in) :: filename
        class(magnetic_field_t), allocatable, intent(out) :: field

        character(len(filename)) :: stripped_name
        type(coils_field_t) :: raw_coils
        type(splined_field_t), allocatable :: splined_coils
        type(vmec_field_t) :: vmec_field
        integer :: file_type, ierr
        character(len=2048) :: message
#ifdef GVEC_AVAILABLE
        class(gvec_field_t), allocatable :: gvec_temp
#endif

        stripped_name = strip_directory(filename)

        if (endswith(filename, '.nc')) then
            call detect_refcoords_file_type(filename, file_type, ierr, message)
            if (ierr /= 0) then
                print *, 'field_from_file: NetCDF file detection error for ', &
                    trim(filename)
                print *, trim(message)
                error stop
            end if

            select case (file_type)
            case (refcoords_file_vmec_wout)
                call create_vmec_field(vmec_field)
                allocate (field, source=vmec_field)
            case (refcoords_file_chartmap)
                print *, &
                    'field_from_file: chartmap NetCDF is a coordinate system file,', &
                    ' not a field:'
                print *, '  filename = ', trim(filename)
                print *, &
                    'Set coord_input to this chartmap file and set field_input to a', &
                    ' VMEC wout.'
                error stop
            case (refcoords_file_unknown)
                print *, 'field_from_file: Unknown NetCDF file type: ', trim(filename)
                error stop
            case default
                print *, 'field_from_file: Unexpected file_type ', file_type, ' for ', &
                    trim(filename)
                error stop
            end select
        else if (startswidth(stripped_name, 'coils') .or. &
                 endswith(filename, '.coils')) then
            call create_coils_field(filename, raw_coils)
            allocate (splined_coils)
            call create_splined_field(raw_coils, ref_coords, splined_coils)
            call move_alloc(splined_coils, field)
        else if (endswith(filename, '.dat')) then
#ifdef GVEC_AVAILABLE
            call create_gvec_field(filename, gvec_temp)
            call move_alloc(gvec_temp, field)
#else
            print *, 'ERROR: GVEC support not compiled. Rebuild with -DENABLE_GVEC=ON'
            error stop
#endif
        else
            print *, 'field_from_file: Unknown file name format ', filename
            error stop
        end if
    end subroutine field_from_file

    function startswidth(text, start)
        logical :: startswidth
        character(*), intent(in) :: text
        character(*), intent(in) :: start
        integer :: len_text, len_start

        len_text = len_trim(text)
        len_start = len_trim(start)

        startswidth = .false.
        if (len_text >= len_start) then
            if (text(1:len_start) == start) then
                startswidth = .true.
            end if
        end if
    end function startswidth

    function endswith(text, ending)
        logical :: endswith
        character(*), intent(in) :: text
        character(*), intent(in) :: ending
        integer :: len_text, len_end

        len_text = len_trim(text)
        len_end = len_trim(ending)

        endswith = .false.
        if (len_text >= len_end) then
            if (text(len_text - len_end + 1:len_text) == ending) then
                endswith = .true.
            end if
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

end module field
