module simple_profiles
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    integer, parameter :: MAX_POWER_SERIES = 10

    integer, parameter :: POWER_SERIES = 1
    integer, parameter :: TWO_POWER = 2

    character(len=32) :: profile_type = "power_series"
    integer :: active_profile = POWER_SERIES

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

    namelist /profiles/ profile_type, &
        Te_scale, Ti1_scale, Ti2_scale, ni1_scale, ni2_scale, &
        Te_coef, Ti1_coef, Ti2_coef, ni1_coef, ni2_coef, &
        Te_p1, Te_p2, Ti1_p1, Ti1_p2, Ti2_p1, Ti2_p2, &
        ni1_p1, ni1_p2, ni2_p1, ni2_p2

contains

    subroutine read_profiles(config_file)
        character(len=*), intent(in) :: config_file
        integer :: ios, unit

        open (newunit=unit, file=config_file, status='old', action='read', iostat=ios)
        if (ios /= 0) return

        read (unit, nml=profiles, iostat=ios)
        close (unit)

        select case (trim(profile_type))
        case ("two_power")
            active_profile = TWO_POWER
        case default
            active_profile = POWER_SERIES
        end select
    end subroutine read_profiles

    subroutine get_plasma_params(s, Te, Ti1, Ti2, ni1, ni2)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: Te, Ti1, Ti2, ni1, ni2

        if (active_profile == TWO_POWER) then
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
