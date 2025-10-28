program test_field_analytical_geoflux

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field, only: field_from_file, MagneticField
    use field_geoflux, only: GeofluxField, geoflux_ready

    implicit none

    real(dp), parameter :: pi_val = acos(-1.0_dp)
    class(MagneticField), allocatable :: field_obj
    real(dp) :: x(3), Acov(3), hcov(3), Bmod

    call field_from_file('analytical', field_obj, .true.)

    if (.not. geoflux_ready) then
        error stop 'Analytical geoflux field failed to initialize'
    end if

    select type(field_obj)
    type is (GeofluxField)
        x = [sqrt(0.25_dp), 0.3_dp, 0.1_dp]
        call field_obj%evaluate(x, Acov, hcov, Bmod)

        if (Bmod <= 0.0_dp) then
            error stop 'Analytical geoflux |B| must be positive'
        end if

        if (abs(Acov(3)) <= 0.0_dp) then
            error stop 'Analytical geoflux Aphi expected non-zero away from axis'
        end if

        if (abs(hcov(3)) <= 1.0e-8_dp) then
            error stop 'Analytical geoflux hphi unexpectedly small'
        end if

        call check_toroidal_component(field_obj)
    class default
        error stop 'field_from_file(analytical) did not return GeofluxField'
    end select

contains

    subroutine check_toroidal_component(field)
        class(GeofluxField), intent(in) :: field
        real(dp) :: min_hp, max_hr
        real(dp) :: r, th, phi
        real(dp) :: Acov_local(3), hcov_local(3), Bmod_local
        integer :: ir, it, ip

        min_hp = huge(1.0_dp)
        max_hr = 0.0_dp

        do ip = 0, 7
            phi = real(ip, dp) * (0.25_dp*pi_val)
            do it = 0, 7
                th = real(it, dp) * (0.25_dp*pi_val)
                do ir = 0, 8
                    r = 1.0e-3_dp + real(ir, dp) * (0.999_dp/8.0_dp)
                    call field%evaluate([r, th, phi], Acov_local, hcov_local, &
                        Bmod_local)
                    min_hp = min(min_hp, abs(hcov_local(3)))
                    max_hr = max(max_hr, abs(hcov_local(1)))
                end do
            end do
        end do

        if (min_hp <= 1.0e-8_dp) then
            error stop 'Analytical geoflux toroidal component vanishes on sample grid'
        end if

        if (max_hr/min_hp > 5.0e-2_dp) then
            error stop 'Analytical geoflux radial component too large relative to toroidal'
        end if
    end subroutine check_toroidal_component

end program test_field_analytical_geoflux
