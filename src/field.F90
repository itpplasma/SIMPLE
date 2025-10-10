module field

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: MagneticField
use field_vmec, only: VmecField
use field_geoflux, only: GeofluxField, initialize_geoflux_field, &
    mark_geoflux_initialized
use analytical_geoflux_field, only: init_analytical_geoflux
use field_coils, only: CoilsField, create_coils_field
#ifdef GVEC_AVAILABLE
use field_gvec, only: GvecField, create_gvec_field
#endif

implicit none

contains

subroutine field_from_file(filename, field)
    use tokamak_config_mod, only: tok_R0, tok_epsilon, tok_kappa, tok_delta, &
        tok_A_param, tok_B0, tok_Nripple, tok_a0, tok_alpha0, tok_delta0, &
        tok_z0
    character(*), intent(in) :: filename
    class(MagneticField), allocatable, intent(out) :: field

    character(len(filename)) :: stripped_name
    character(len(filename)) :: lower_name
    logical :: use_analytical_field
    class(CoilsField), allocatable :: coils_temp
#ifdef GVEC_AVAILABLE
    class(GvecField), allocatable :: gvec_temp
#endif

    stripped_name = strip_directory(filename)
    lower_name = to_lower(filename)
    use_analytical_field = index(lower_name, 'analytical') > 0 .or. &
        index(lower_name, 'tokamak') > 0

    if (is_geqdsk(filename)) then
        call initialize_geoflux_field(trim(filename))
        allocate(GeofluxField :: field)
    else if (use_analytical_field) then
        call init_analytical_geoflux(tok_R0, tok_epsilon, tok_kappa, tok_delta, &
            tok_A_param, tok_B0, tok_Nripple, tok_a0, tok_alpha0, tok_delta0, &
            tok_z0)
        call mark_geoflux_initialized(trim(lower_name), .true.)
        allocate(GeofluxField :: field)
    else if (endswith(filename, '.nc')) then
        allocate(VmecField :: field)
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


logical function is_geqdsk(filename)
    character(*), intent(in) :: filename

    character(:), allocatable :: lower_name

    lower_name = to_lower(trim(filename))

    is_geqdsk = endswith(lower_name, '.geqdsk') .or. endswith(lower_name, '.eqdsk')
    if (.not. is_geqdsk) then
        is_geqdsk = startswidth(strip_directory(lower_name), 'geqdsk')
    end if
end function is_geqdsk


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


end module field
