program test_spectre_fo_spline
    ! Independent raw-field oracle at off-grid points, including both traces
    ! of every sheet. Grid refinement must improve the interpolated field.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_spectre, only: spectre_field_t, create_spectre_field
    use magfie_sub, only: set_magfie_spectre_field, spectre_field, &
        TESLA_TO_GAUSS, M_TO_CM
    use spectre_fo_field, only: spectre_fo_field_build, spectre_fo_field_local, &
        spectre_fo_field_locate, spectre_fo_local_t
    implicit none
    type(spectre_field_t) :: field
    type(spectre_fo_local_t) :: loc
    character(len=1024) :: path
    real(dp) :: u(3), raw(3), a(3), hc(3), bmag, sqgb(3), e(3, 3), g(3, 3)
    real(dp) :: gi(3, 3), jac, b(3), located(3), errors(3, 2), fraction, delta
    integer :: status, level, volume, j, owner, nr, nt, nz

    call get_command_argument(1, path)
    call create_spectre_field(field, trim(path), status)
    if (status /= 0) error stop 'spline field initialization'
    call set_magfie_spectre_field(field)
    errors = 0.0_dp
    do level = 1, 2
        nr = 16*level
        nt = 32*level
        nz = 16*level
        if (all(field%data%in == 0)) nz = 7
        call spectre_fo_field_build(nr, nt, nz)
        do volume = 1, field%data%Mvol
            do j = 1, 48
                fraction = 0.1_dp + 0.8_dp*real(mod(13*j, 47), dp)/47
                if (j == 1) fraction = 1.0e-10_dp
                if (j == 2) fraction = 1.0_dp - 1.0e-10_dp
                if (volume == 1 .and. j == 1) fraction = 0.05_dp
                u = [real(volume - 1, dp) + fraction, 0.713_dp*j, 0.317_dp*j]
                call spectre_fo_field_local(u, volume, loc)
                call spectre_field%coords%evaluate_cart(u, raw)
                errors(1, level) = max(errors(1, level), &
                    norm2(loc%x/M_TO_CM - raw)/norm2(raw))
                call spectre_field%evaluate(u, a, hc, bmag, sqgb)
                call spectre_field%coords%covariant_basis(u, e)
                call spectre_field%coords%metric_tensor(u, g, gi, jac)
                b = matmul(e, sqgb)/jac
                errors(2, level) = max(errors(2, level), &
                    norm2(loc%b*loc%bmod/TESLA_TO_GAUSS - b)/bmag)
                errors(3, level) = max(errors(3, level), &
                    norm2(loc%acov/(TESLA_TO_GAUSS*M_TO_CM**2) - a) &
                    /max(norm2(a), 1.0e-6_dp))
                located = u + [0.001_dp, 0.002_dp, 0.0_dp]
                owner = volume
                call spectre_fo_field_locate(loc%x, located, owner, status)
                delta = norm2(located - u)
                if (status /= 0 .or. owner /= volume .or. delta > 1.0e-8_dp) then
                    print *, 'inverse mismatch', level, volume, j, status, owner
                    print *, 'target', u, 'located', located, 'delta', delta
                    error stop 'spline coordinate inverse'
                end if
            end do
        end do
        print '(a,i1,3es14.5)', 'grid level ', level, errors(:, level)
    end do
    if (errors(1, 2) > 1.0e-7_dp) error stop 'spline geometry error'
    if (errors(2, 2) > 1.0e-5_dp) error stop 'spline magnetic field error'
    if (errors(3, 2) > 1.0e-6_dp) error stop 'spline potential error'
    if (errors(2, 1) > 1.0e-10_dp) then
        if (errors(2, 2) > 0.25_dp*errors(2, 1)) &
            error stop 'field did not converge with spline refinement'
    end if
end program test_spectre_fo_spline
