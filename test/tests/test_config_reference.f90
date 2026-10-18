program test_config_reference
    !> The documented input reference must name every &config and &profiles
    !> value with its compiled default, reset_config must restore those
    !> defaults, and input files using removed keys must still be read.
    use params, only: config, reset_config, read_config_namelist, read_config, &
        trace_time, ntestpart
    use simple_profiles, only: profiles, reset_profiles, read_profiles, Te_scale

    implicit none

    character(*), parameter :: removed_keys_file = 'removed_keys.in'
    character(1024) :: reference_file
    character(1024), allocatable :: defaults(:), documented(:)
    integer :: nerr = 0

    if (command_argument_count() /= 1) error stop 'usage: test_config_reference.x FILE'
    call get_command_argument(1, reference_file)

    call reset_config
    call reset_profiles
    call dump_namelists(defaults)

    call read_config_namelist(trim(reference_file))
    call read_profiles(trim(reference_file))
    call dump_namelists(documented)
    call compare(defaults, documented, 'reference file vs compiled defaults')
    call check_names_written(defaults, trim(reference_file))

    trace_time = 42.0d0
    Te_scale = 42.0d0
    call reset_config
    call reset_profiles
    call dump_namelists(documented)
    call compare(defaults, documented, 'reset vs compiled defaults')

    call write_removed_keys_file
    call read_config(removed_keys_file)
    if (ntestpart /= 3) call fail('file with removed keys was not read')

    if (nerr > 0) error stop 'test_config_reference failed'
    print *, 'test_config_reference passed'

contains

    subroutine dump_namelists(lines)
        character(1024), allocatable, intent(out) :: lines(:)
        character(1024) :: buffer(4096)
        integer :: unit, n, ios

        open (newunit=unit, status='scratch', delim='apostrophe')
        write (unit, nml=config)
        write (unit, nml=profiles)
        rewind (unit)
        n = 0
        do
            read (unit, '(a)', iostat=ios) buffer(n + 1)
            if (ios /= 0) exit
            n = n + 1
        end do
        close (unit)
        lines = buffer(1:n)
    end subroutine dump_namelists

    subroutine compare(expected, actual, what)
        character(*), intent(in) :: expected(:), actual(:), what
        integer :: i

        if (size(expected) /= size(actual)) then
            call fail(what//': different number of namelist lines')
            return
        end if
        do i = 1, size(expected)
            if (expected(i) /= actual(i)) then
                call fail(what//':')
                print *, '  expected: ', trim(expected(i))
                print *, '  found:    ', trim(actual(i))
            end if
        end do
    end subroutine compare

    subroutine check_names_written(dump, path)
        ! Each name must be assigned explicitly so that it is documented.
        character(*), intent(in) :: dump(:), path
        character(1024) :: line, name
        logical :: found
        integer :: i, eq, unit, ios

        do i = 1, size(dump)
            eq = index(dump(i), '=')
            if (eq == 0 .or. dump(i) (1:1) /= ' ') cycle
            name = adjustl(to_lower(dump(i) (1:eq - 1)))
            found = .false.
            open (newunit=unit, file=path, status='old', action='read')
            do
                read (unit, '(a)', iostat=ios) line
                if (ios /= 0) exit
                line = adjustl(to_lower(line))
                if (index(line, trim(name)//' ') == 1 .or. &
                    index(line, trim(name)//'=') == 1 .or. &
                    index(line, trim(name)//'(') == 1) then
                    found = .true.
                    exit
                end if
            end do
            close (unit)
            if (.not. found) call fail('not documented: '//trim(name))
        end do
    end subroutine check_names_written

    subroutine write_removed_keys_file
        integer :: unit

        open (newunit=unit, file=removed_keys_file, status='replace')
        write (unit, '(a)') '&config', 'ntestpart = 3', 'local = .false.', &
            'loopskip = 2', 'bmod_ref = 3d4', 'macrostep_log_ratio = 10d0', '/'
        close (unit)
    end subroutine write_removed_keys_file

    pure function to_lower(text) result(lower)
        character(*), intent(in) :: text
        character(len(text)) :: lower
        integer :: i, code

        lower = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) &
                lower(i:i) = achar(code + 32)
        end do
    end function to_lower

    subroutine fail(message)
        character(*), intent(in) :: message

        print *, 'FAIL: ', message
        nerr = nerr + 1
    end subroutine fail

end program test_config_reference
