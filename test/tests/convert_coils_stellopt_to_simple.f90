program convert_coils_stellopt_to_simple
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit

    implicit none

    character(len=2048) :: input_file, output_file, line
    integer :: in_unit, out_unit, ios, n, capacity, i
    logical :: in_filament

    real(dp), allocatable :: x(:), y(:), z(:), cur(:)
    real(dp) :: xv, yv, zv, cv

    call get_command_argument(1, input_file)
    call get_command_argument(2, output_file)

    if (len_trim(input_file) == 0 .or. len_trim(output_file) == 0) then
        write(error_unit, '(a)') 'Usage: convert_coils_stellopt_to_simple.x ' // &
                                 '<input.coils> <output.simple>'
        error stop 1
    end if

    capacity = 4096
    allocate(x(capacity), y(capacity), z(capacity), cur(capacity))
    n = 0
    in_filament = .false.

    open(newunit=in_unit, file=trim(input_file), status='old', action='read')
    do
        read(in_unit, '(a)', iostat=ios) line
        if (ios /= 0) exit

        call to_lower_inplace(line)
        if (index(line, 'begin filament') > 0) then
            in_filament = .true.
            cycle
        end if

        if (.not. in_filament) cycle
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle

        read(line, *, iostat=ios) xv, yv, zv, cv
        if (ios /= 0) cycle

        if (n == capacity) then
            call grow_arrays(x, y, z, cur, capacity)
        end if

        n = n + 1
        x(n) = xv
        y(n) = yv
        z(n) = zv
        cur(n) = cv
    end do
    close(in_unit)

    open(newunit=out_unit, file=trim(output_file), status='replace', action='write')
    write(out_unit, '(i0)') n
    do i = 1, n
        write(out_unit, '(4(es24.16,1x))') x(i), y(i), z(i), cur(i)
    end do
    close(out_unit)

contains

    subroutine to_lower_inplace(s)
        character(len=*), intent(inout) :: s
        integer :: i, c

        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) then
                s(i:i) = achar(c + (iachar('a') - iachar('A')))
            end if
        end do
    end subroutine to_lower_inplace

    subroutine grow_arrays(x, y, z, cur, capacity)
        real(dp), allocatable, intent(inout) :: x(:), y(:), z(:), cur(:)
        integer, intent(inout) :: capacity

        real(dp), allocatable :: xn(:), yn(:), zn(:), cn(:)

        capacity = capacity * 2
        allocate(xn(capacity), yn(capacity), zn(capacity), cn(capacity))

        xn(1:size(x)) = x
        yn(1:size(y)) = y
        zn(1:size(z)) = z
        cn(1:size(cur)) = cur

        call move_alloc(xn, x)
        call move_alloc(yn, y)
        call move_alloc(zn, z)
        call move_alloc(cn, cur)
    end subroutine grow_arrays

end program convert_coils_stellopt_to_simple
