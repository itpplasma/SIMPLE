module field

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: MagneticField
use field_vmec, only: VmecField
use field_coils, only: CoilsField, create_coils_field

implicit none

contains

function field_from_file(filename)
    class(MagneticField), allocatable :: field_from_file
    character(*), intent(in) :: filename

    character(len(filename)) :: stripped_name
    stripped_name = strip_directory(filename)

    if (endswith(filename, '.nc')) then
        allocate(VmecField :: field_from_file)
    else if (startswidth(stripped_name, 'coils') .or. endswith(filename, '.coils')) then
        field_from_file = create_coils_field(filename)
    else
        print *,  'field_from_file: Unknown file name format ', filename
        error stop
    end if
end function field_from_file


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
