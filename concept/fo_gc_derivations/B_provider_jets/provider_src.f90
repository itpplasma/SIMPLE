module provider_src
    !> Read a right-hand side from a cited physical line of a Fortran source,
    !> normalise its literals and rename cited array references, then parse it.
    !> Only the operator structure of the actual source is used; the renames are
    !> explicit (array element text -> symbol) and recorded by the caller.
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t
    use fortsym_dialect, only: dialect, DIA_FORTRAN
    use fortsym_parse, only: parse_expr_in
    implicit none
    private
    public :: source_rhs, parse_rhs, replace_all, normalise_literals, strip_real_dp
    public :: source_line

contains

    !> Physical line `line` of `file` (1-based), without trailing blanks.
    function source_line(file, line) result(text)
        character(*), intent(in) :: file
        integer, intent(in) :: line
        character(:), allocatable :: text
        character(4096) :: raw
        integer :: u, ios, k
        open (newunit=u, file=file, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'provider_src: cannot open source'
        do k = 1, line
            read (u, '(a)', iostat=ios) raw
            if (ios /= 0) error stop 'provider_src: line beyond end of file'
        end do
        close (u)
        text = trim(raw)
    end function source_line

    !> Right-hand side of the assignment starting at physical line `line`,
    !> continuation lines joined and the trailing comment removed.
    function source_rhs(file, line) result(rhs)
        character(*), intent(in) :: file
        integer, intent(in) :: line
        character(:), allocatable :: rhs, text, piece
        integer :: k, eq, depth, i
        text = ''
        k = line
        do
            piece = strip_comment(source_line(file, k))
            piece = trim(adjustl(piece))
            if (len(text) > 0 .and. len(piece) > 0) then
                if (piece(1:1) == '&') piece = piece(2:)
            end if
            if (len(piece) > 0) then
                if (piece(len(piece):len(piece)) == '&') then
                    text = text//' '//piece(1:len(piece) - 1)
                    k = k + 1
                    cycle
                end if
            end if
            text = text//' '//piece
            exit
        end do
        eq = 0
        depth = 0
        do i = 1, len(text)
            select case (text(i:i))
            case ('(', '[')
                depth = depth + 1
            case (')', ']')
                depth = depth - 1
            case ('=')
                if (depth == 0 .and. eq == 0) then
                    if (i < len(text)) then
                        if (text(i + 1:i + 1) == '=') cycle
                    end if
                    if (i > 1) then
                        if (index('/<>=', text(i - 1:i - 1)) > 0) cycle
                    end if
                    eq = i
                end if
            end select
        end do
        if (eq == 0) error stop 'provider_src: no assignment on cited line'
        rhs = trim(adjustl(text(eq + 1:)))
    end function source_rhs

    function strip_comment(text) result(r)
        character(*), intent(in) :: text
        character(:), allocatable :: r
        integer :: i
        r = text
        i = index(text, '!')
        if (i > 0) r = text(1:i - 1)
    end function strip_comment

    function replace_all(text, old, new) result(r)
        character(*), intent(in) :: text, old, new
        character(:), allocatable :: r
        integer :: i, start
        r = ''
        start = 1
        do
            i = index(text(start:), old)
            if (i == 0) exit
            r = r//text(start:start + i - 2)//new
            start = start + i - 1 + len(old)
            if (start > len(text)) exit
        end do
        if (start <= len(text)) r = r//text(start:)
    end function replace_all

    !> real(X, dp) -> (X); exact conversion of an integer expression.
    function strip_real_dp(text) result(r)
        character(*), intent(in) :: text
        character(:), allocatable :: r
        integer :: i, j, depth, comma
        r = text
        do
            i = index(r, 'real(')
            if (i == 0) exit
            if (i > 1) then
                if (verify(r(i - 1:i - 1), &
                    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_') == 0) &
                    error stop 'provider_src: real( inside identifier'
            end if
            depth = 0
            comma = 0
            do j = i + 4, len(r)
                if (r(j:j) == '(') depth = depth + 1
                if (r(j:j) == ')') depth = depth - 1
                if (r(j:j) == ',' .and. depth == 1) comma = j
                if (depth == 0) exit
            end do
            if (comma == 0) error stop 'provider_src: real() without kind'
            if (trim(adjustl(r(comma + 1:j - 1))) /= 'dp') error stop 'provider_src: kind'
            r = r(1:i - 1)//'('//r(i + 5:comma - 1)//')'//r(j + 1:)
        end do
    end function strip_real_dp

    !> Decimal literals d.d_dp -> exact rational (N/D). Every literal in the
    !> cited lines is a short decimal; the conversion is exact for those that
    !> are exactly representable (the caller only cites dyadic literals).
    function normalise_literals(text) result(r)
        character(*), intent(in) :: text
        character(:), allocatable :: r, ip, fp
        character(64) :: buf
        integer :: i, j, k, p, e10
        integer(8) :: num, den
        r = ''
        i = 1
        do while (i <= len(text))
            k = index(text(i:), '_dp')
            if (k == 0) then
                r = r//text(i:)
                exit
            end if
            k = i + k - 1
            p = k - 1
            do while (p >= i)
                if (index('0123456789.', text(p:p)) == 0) exit
                p = p - 1
            end do
            j = index(text(p + 1:k - 1), '.')
            if (j == 0) error stop 'provider_src: literal without point'
            ip = text(p + 1:p + j - 1)
            fp = text(p + j + 1:k - 1)
            if (len(ip) == 0) ip = '0'
            e10 = len(fp)
            den = 10_8**e10
            buf = ip//fp
            read (buf, *) num
            write (buf, '(a,i0,a,i0,a)') '(', num, '/', den, ')'
            r = r//text(i:p)//trim(buf)
            i = k + 3
        end do
    end function normalise_literals

    function parse_rhs(a, rhs, olds, news) result(e)
        type(arena_t), target, intent(inout) :: a
        character(*), intent(in) :: rhs
        character(*), intent(in) :: olds(:), news(:)
        type(expr_t) :: e
        character(:), allocatable :: t, msg
        logical :: ok
        integer :: k
        t = rhs
        do k = 1, size(olds)
            t = replace_all(t, trim(olds(k)), trim(news(k)))
        end do
        t = strip_real_dp(t)
        t = normalise_literals(t)
        e = parse_expr_in(a, t, dialect(DIA_FORTRAN), ok, msg)
        if (.not. ok) then
            print *, 'parse refused: ', t
            print *, msg
            error stop 'provider_src: parse'
        end if
    end function parse_rhs
end module provider_src
