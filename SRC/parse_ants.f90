module parse_ants
implicit none
integer, parameter :: maxlen = 350

contains

subroutine process_line(line, v_par, v_perp, u, v, s)
    character(maxlen), intent(in) :: line
    integer :: i
    real(8), intent(out) :: v_par, v_perp, u, v, s

    do i = 1, maxlen
        if (line(i-6:i) == 'v_par =') read(line(i+1:maxlen), *) v_par
        if (line(i-7:i) == 'v_perp =') read(line(i+1:maxlen), *) v_perp
        if (line(i-2:i) == 'u =') read(line(i+1:maxlen), *) u
        if (line(i-2:i) == 'v =') read(line(i+1:maxlen), *) v
        if (line(i-2:i) == 's =') read(line(i+1:maxlen), *) s
    end do
end subroutine process_line
end module parse_ants
