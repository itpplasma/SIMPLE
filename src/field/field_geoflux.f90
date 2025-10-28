module field_geoflux

use, intrinsic :: iso_fortran_env, only: dp => real64
use field_base, only: MagneticField
use geoflux_field, only: spline_geoflux_data, splint_geoflux_field
use analytical_geoflux_field, only: splint_analytical_geoflux_field

implicit none

type, extends(MagneticField) :: GeofluxField
contains
    procedure :: evaluate => geoflux_evaluate
end type GeofluxField

character(len=:), allocatable :: cached_geqdsk
logical :: geoflux_ready = .false.
integer, parameter :: default_ns_cache = 128
integer, parameter :: default_ntheta_cache = 256
logical :: analytical_mode = .false.

contains

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

    call mark_geoflux_initialized(normalized_path, .false.)
end subroutine initialize_geoflux_field


subroutine mark_geoflux_initialized(cache_key, use_analytical)
    character(len=*), intent(in), optional :: cache_key
    logical, intent(in), optional :: use_analytical

    if (present(cache_key)) then
        cached_geqdsk = trim(adjustl(cache_key))
    end if

    analytical_mode = .false.
    if (present(use_analytical)) analytical_mode = use_analytical

    geoflux_ready = .true.
end subroutine mark_geoflux_initialized


subroutine geoflux_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    class(GeofluxField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: Acov(3)
    real(dp), intent(out) :: hcov(3)
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    real(dp) :: geo_s, theta_geo, phi_geo
    real(dp) :: Acov_geo(3), hcov_geo(3)
    real(dp) :: sqgBctr_geo(3)
    real(dp) :: ds_dr

    if (.not. geoflux_ready) then
        error stop 'GeofluxField: call initialize_geoflux_field before evaluation'
    end if

    geo_s = max(0.0_dp, min(1.0_dp, x(1)*x(1)))
    theta_geo = x(2)
    phi_geo = x(3)
    ds_dr = 2.0_dp * x(1)

    if (analytical_mode) then
        if (present(sqgBctr)) then
            call splint_analytical_geoflux_field(geo_s, theta_geo, phi_geo, &
                Acov_geo, hcov_geo, Bmod, sqgBctr_geo)
        else
            call splint_analytical_geoflux_field(geo_s, theta_geo, phi_geo, &
                Acov_geo, hcov_geo, Bmod)
        end if
    else
        if (present(sqgBctr)) then
            call splint_geoflux_field(geo_s, theta_geo, phi_geo, Acov_geo, hcov_geo, &
                Bmod, sqgBctr_geo)
        else
            call splint_geoflux_field(geo_s, theta_geo, phi_geo, Acov_geo, hcov_geo, &
                Bmod)
        end if
    end if

    Acov = Acov_geo
    hcov = hcov_geo

    Acov(1) = Acov(1) * ds_dr
    hcov(1) = hcov(1) * ds_dr

    if (present(sqgBctr)) then
        sqgBctr = sqgBctr_geo
        sqgBctr(1) = sqgBctr(1) * ds_dr
    end if
end subroutine geoflux_evaluate


logical function geoflux_is_analytical()
    geoflux_is_analytical = analytical_mode
end function geoflux_is_analytical

end module field_geoflux
