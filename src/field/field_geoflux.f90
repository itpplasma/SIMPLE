module field_geoflux

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: magnetic_field_t
    use libneo_coordinates, only: make_geoflux_coordinate_system
    use geoflux_field, only: spline_geoflux_data, splint_geoflux_field

    implicit none

    type, extends(magnetic_field_t) :: geoflux_field_t
    contains
        procedure :: evaluate => geoflux_evaluate
    end type geoflux_field_t

    character(len=:), allocatable :: cached_geqdsk
    logical :: geoflux_ready = .false.
    integer, parameter :: default_ns_cache = 128
    integer, parameter :: default_ntheta_cache = 256

contains

    subroutine create_geoflux_field(field)
        type(geoflux_field_t), intent(out) :: field

        call make_geoflux_coordinate_system(field%coords)
    end subroutine create_geoflux_field

    subroutine initialize_geoflux_field(geqdsk_file, ns_cache, ntheta_cache)
        character(len=*), intent(in) :: geqdsk_file
        integer, intent(in), optional :: ns_cache, ntheta_cache
        integer :: ns_val, ntheta_val
        character(len=:), allocatable :: normalized_path

        ns_val = default_ns_cache
        if (present(ns_cache)) ns_val = ns_cache

        ntheta_val = default_ntheta_cache
        if (present(ntheta_cache)) ntheta_val = ntheta_cache

        normalized_path = trim(adjustl(geqdsk_file))

        if (geoflux_ready) then
            if (allocated(cached_geqdsk)) then
                if (normalized_path == cached_geqdsk) return
            end if
        end if

        call spline_geoflux_data(normalized_path, ns_val, ntheta_val)

        cached_geqdsk = normalized_path
        geoflux_ready = .true.
    end subroutine initialize_geoflux_field

    subroutine geoflux_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        class(geoflux_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3)
        real(dp), intent(out) :: hcov(3)
        real(dp), intent(out) :: Bmod
        real(dp), intent(out), optional :: sqgBctr(3)

        real(dp) :: geo_s, theta_geo, phi_geo
        real(dp) :: Acov_geo(3), hcov_geo(3)
        real(dp) :: sqgBctr_geo(3)

        if (.not. geoflux_ready) then
            error stop &
                'geoflux_field_t: call initialize_geoflux_field before evaluation'
        end if

        geo_s = max(0.0_dp, min(1.0_dp, x(1)))
        theta_geo = x(2)
        phi_geo = x(3)

        if (present(sqgBctr)) then
            call splint_geoflux_field(geo_s, theta_geo, phi_geo, Acov_geo, &
                hcov_geo, Bmod, &
                                      sqgBctr_geo)
        else
            call splint_geoflux_field(geo_s, theta_geo, phi_geo, Acov_geo, &
                hcov_geo, Bmod)
        end if

        Acov = Acov_geo
        hcov = hcov_geo

        if (present(sqgBctr)) then
            sqgBctr = sqgBctr_geo
        end if
    end subroutine geoflux_evaluate

end module field_geoflux
