module field

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: magnetic_field_t
use field_vmec, only: vmec_field_t
use field_coils, only: coils_field_t, create_coils_field
#ifdef GVEC_AVAILABLE
use field_gvec, only: gvec_field_t, create_gvec_field
#endif

implicit none

contains

subroutine field_from_file(filename, field)
    character(*), intent(in) :: filename
    class(magnetic_field_t), allocatable, intent(out) :: field

    character(len(filename)) :: stripped_name
    class(coils_field_t), allocatable :: coils_temp
#ifdef GVEC_AVAILABLE
    class(gvec_field_t), allocatable :: gvec_temp
#endif

    stripped_name = strip_directory(filename)

    if (endswith(filename, '.nc')) then
        allocate(vmec_field_t :: field)
    else if (startswidth(stripped_name, 'coils') .or. endswith(filename, '.coils')) then
        call create_coils_field(filename, coils_temp)
        call move_alloc(coils_temp, field)
    else if (endswith(filename, '.dat')) then
#ifdef GVEC_AVAILABLE
        call create_gvec_field(filename, gvec_temp)
        call move_alloc(gvec_temp, field)
#else
        print *, 'ERROR: GVEC support not compiled. Rebuild with -DENABLE_GVEC=ON'
        error stop
#endif
    else
        print *,  'field_from_file: Unknown file name format ', filename
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
        if (text(1 : len_start) == start) then
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
        if (text(len_text - len_end + 1 : len_text) == ending) then
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
            strip_directory = filename(i + 1 : len(filename))
            return
        end if
    end do
end function strip_directory


end module field
