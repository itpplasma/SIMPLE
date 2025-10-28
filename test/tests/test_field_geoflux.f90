program test_field_geoflux

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field, only: field_from_file, MagneticField
    use field_geoflux, only: GeofluxField

    implicit none

    class(MagneticField), allocatable :: field_obj
    real(dp) :: x(3), Acov(3), hcov(3), Bmod
    real(dp) :: norm_sq
    character(len=512) :: geqdsk_path
    integer :: status

    geqdsk_path = 'EQDSK_I.geqdsk'
    call get_environment_variable('LIBNEO_TEST_GEQDSK', value=geqdsk_path, status=status)
    if (status /= 0 .or. len_trim(geqdsk_path) == 0) then
        geqdsk_path = 'EQDSK_I.geqdsk'
    end if

    call field_from_file(trim(geqdsk_path), field_obj)

    select type(field_obj)
    type is (GeofluxField)
        x = [sqrt(0.25_dp), 0.3_dp, 0.0_dp]
        call field_obj%evaluate(x, Acov, hcov, Bmod)

        if (Bmod <= 0.0_dp) then
            error stop 'GeofluxField: Bmod must be positive'
        end if

        if (abs(Acov(3)) <= 0.0_dp) then
            error stop 'GeofluxField: Aphi expected to be non-zero away from axis'
        end if

        norm_sq = sum(hcov**2)
        if (norm_sq <= 0.0_dp) then
            error stop 'GeofluxField: hcov has zero magnitude'
        end if
    class default
        error stop 'field_from_file did not return GeofluxField for GEQDSK input'
    end select

end program test_field_geoflux
