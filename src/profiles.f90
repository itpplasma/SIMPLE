module simple_profiles
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    integer, parameter :: MAX_POWER_SERIES = 10

    integer, parameter :: POWER_SERIES = 1
    integer, parameter :: TWO_POWER = 2

    character(len=32) :: profile_type = "power_series"
    character(len=1000) :: profile_file = ""
    integer :: active_profile = POWER_SERIES
    real(dp), allocatable, private :: table_radius(:), table_values(:, :)

    real(dp) :: Te_scale = 1.0d4
    real(dp) :: Ti1_scale = 1.0d4
    real(dp) :: Ti2_scale = 1.0d4
    real(dp) :: ni1_scale = 0.5d20
    real(dp) :: ni2_scale = 0.5d20

    real(dp) :: Te_coef(MAX_POWER_SERIES) = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    real(dp) :: Ti1_coef(MAX_POWER_SERIES) = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    real(dp) :: Ti2_coef(MAX_POWER_SERIES) = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    real(dp) :: ni1_coef(MAX_POWER_SERIES) = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    real(dp) :: ni2_coef(MAX_POWER_SERIES) = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]

    real(dp) :: Te_p1 = 1.0d0, Te_p2 = 0.0d0
    real(dp) :: Ti1_p1 = 1.0d0, Ti1_p2 = 0.0d0
    real(dp) :: Ti2_p1 = 1.0d0, Ti2_p2 = 0.0d0
    real(dp) :: ni1_p1 = 1.0d0, ni1_p2 = 0.0d0
    real(dp) :: ni2_p1 = 1.0d0, ni2_p2 = 0.0d0

    namelist /profiles/ profile_type, profile_file, &
        Te_scale, Ti1_scale, Ti2_scale, ni1_scale, ni2_scale, &
        Te_coef, Ti1_coef, Ti2_coef, ni1_coef, ni2_coef, &
        Te_p1, Te_p2, Ti1_p1, Ti1_p2, Ti2_p1, Ti2_p2, &
        ni1_p1, ni1_p2, ni2_p1, ni2_p2

contains

    subroutine export_profile_configuration(configuration, table)
        character(:), allocatable, intent(out) :: configuration
        real(dp), allocatable, intent(out) :: table(:, :)
        character(4096) :: record
        integer :: unit, status

        open (newunit=unit, status='scratch', form='formatted')
        write (unit, nml=profiles)
        rewind (unit)
        configuration = ''
        do
            read (unit, '(a)', iostat=status) record
            if (status < 0) exit
            if (status /= 0) error stop 'Cannot export profile configuration'
            configuration = configuration//trim(record)//new_line('a')
        end do
        close (unit)
        allocate (table(6, 0))
        if (allocated(table_radius)) then
            deallocate (table)
            allocate (table(6, size(table_radius)))
            table(1, :) = table_radius
            table(2:6, :) = table_values
        end if
    end subroutine export_profile_configuration

    subroutine read_profiles(config_file)
        character(len=*), intent(in) :: config_file
        integer :: ios, unit

        open (newunit=unit, file=config_file, status='old', action='read', iostat=ios)
        if (ios /= 0) return

        read (unit, nml=profiles, iostat=ios)
        close (unit)
        if (ios > 0) error stop 'Invalid profiles namelist'
        call select_active_profile
        if (allocated(table_radius)) deallocate (table_radius, table_values)
        if (len_trim(profile_file) > 0) call read_profile_table(profile_file)
    end subroutine read_profiles

    subroutine reset_profiles
        ! Restore every &profiles value to its compiled default. The first call
        ! records the defaults, so call it before changing any value.
        integer, save :: defaults_unit
        logical, save :: recorded = .false.

        if (allocated(table_radius)) deallocate (table_radius, table_values)

        if (.not. recorded) then
            open (newunit=defaults_unit, status='scratch', delim='apostrophe')
            write (defaults_unit, nml=profiles)
            recorded = .true.
        end if
        rewind (defaults_unit)
        read (defaults_unit, nml=profiles)
        call select_active_profile
    end subroutine reset_profiles

    subroutine select_active_profile
        select case (trim(profile_type))
        case ("two_power")
            active_profile = TWO_POWER
        case default
            active_profile = POWER_SERIES
        end select
    end subroutine select_active_profile

    !> Library plasma table: toroidal flux s (0 to smax >= 1), Te, Ti1, Ti2 [eV], ni1, ni2 [m^-3].
    subroutine set_profile_table(n, table, status)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        integer, intent(in) :: n
        real(dp), intent(in) :: table(6, n)
        integer, intent(out) :: status

        status = 1
        if (n < 2) return
        if (.not. all(ieee_is_finite(table))) return
        if (any(table(2:6, :) < 0.0_dp)) return
        if (table(1, 1) /= 0.0_dp .or. table(1, n) < 1.0_dp) return
        if (any(table(1, 2:n) <= table(1, 1:n - 1))) return
        if (allocated(table_radius)) deallocate (table_radius, table_values)
        allocate (table_radius(n), table_values(5, n))
        table_radius = table(1, :)
        table_values = table(2:6, :)
        status = 0
    end subroutine set_profile_table

    subroutine read_profile_table(filename)
        character(len=*), intent(in) :: filename
        real(dp), allocatable :: table(:, :)
        integer :: unit, n, i, ios, status

        open (newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'Cannot open profile_file'
        read (unit, *, iostat=ios) n
        if (ios /= 0) error stop 'Cannot read profile_file row count'
        if (n < 2) error stop 'profile_file needs at least two rows'
        allocate (table(6, n))
        do i = 1, n
            read (unit, *, iostat=ios) table(:, i)
            if (ios /= 0) error stop 'Cannot read profile_file row'
        end do
        close (unit)
        call set_profile_table(n, table, status)
        if (status /= 0) error stop 'Invalid profile_file values or radius grid'
    end subroutine read_profile_table

    function profile_max_s() result(smax)
        real(dp) :: smax
        smax = 1.0_dp
        if (allocated(table_radius)) smax = table_radius(size(table_radius))
    end function profile_max_s

    subroutine get_plasma_params(s, Te, Ti1, Ti2, ni1, ni2)
        use plag_coeff_sub, only: plag_coeff
        real(dp), intent(in) :: s
        real(dp), intent(out) :: Te, Ti1, Ti2, ni1, ni2
        real(dp) :: radius, coefficient(0:0, 2), values(5)
        integer :: i

        if (allocated(table_radius)) then
            radius = max(0.0_dp, min(table_radius(size(table_radius)), s))
            i = max(1, min(size(table_radius) - 1, count(table_radius <= radius)))
            call plag_coeff(2, 0, radius, table_radius(i:i + 1), coefficient)
            values = matmul(table_values(:, i:i + 1), coefficient(0, :))
            Te = values(1)
            Ti1 = values(2)
            Ti2 = values(3)
            ni1 = values(4)
            ni2 = values(5)
        else if (active_profile == TWO_POWER) then
            call get_two_power(s, Te, Ti1, Ti2, ni1, ni2)
        else
            call get_power_series(s, Te, Ti1, Ti2, ni1, ni2)
        end if
    end subroutine get_plasma_params

    pure function eval_power_series(s, scale, coef) result(f)
        real(dp), intent(in) :: s
        real(dp), intent(in) :: scale
        real(dp), intent(in) :: coef(MAX_POWER_SERIES)
        real(dp) :: f
        real(dp) :: s_power
        integer :: n

        f = 0.0d0
        s_power = 1.0d0
        do n = 1, MAX_POWER_SERIES
            f = f + coef(n)*s_power
            s_power = s_power*s
        end do
        f = scale*f
    end function eval_power_series

    pure function eval_two_power(s, scale, p1, p2) result(f)
        real(dp), intent(in) :: s
        real(dp), intent(in) :: scale
        real(dp), intent(in) :: p1
        real(dp), intent(in) :: p2
        real(dp) :: f
        real(dp) :: base

        base = 1.0d0 - s**p1
        if (base < 0.0d0) then
            f = 0.0d0
        else if (abs(base) < epsilon(1.0d0)) then
            if (abs(p2) < epsilon(1.0d0)) then
                f = scale
            else
                f = 0.0d0
            end if
        else
            f = scale*base**p2
        end if
    end function eval_two_power

    subroutine get_power_series(s, Te, Ti1, Ti2, ni1, ni2)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: Te, Ti1, Ti2, ni1, ni2

        Te = eval_power_series(s, Te_scale, Te_coef)
        Ti1 = eval_power_series(s, Ti1_scale, Ti1_coef)
        Ti2 = eval_power_series(s, Ti2_scale, Ti2_coef)
        ni1 = eval_power_series(s, ni1_scale, ni1_coef)
        ni2 = eval_power_series(s, ni2_scale, ni2_coef)
    end subroutine get_power_series

    subroutine get_two_power(s, Te, Ti1, Ti2, ni1, ni2)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: Te, Ti1, Ti2, ni1, ni2

        Te = eval_two_power(s, Te_scale, Te_p1, Te_p2)
        Ti1 = eval_two_power(s, Ti1_scale, Ti1_p1, Ti1_p2)
        Ti2 = eval_two_power(s, Ti2_scale, Ti2_p1, Ti2_p2)
        ni1 = eval_two_power(s, ni1_scale, ni1_p1, ni1_p2)
        ni2 = eval_two_power(s, ni2_scale, ni2_p1, ni2_p2)
    end subroutine get_two_power

end module simple_profiles
