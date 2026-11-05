program test_canflux_axis_continuation
    !> Canonical flux coordinates must use the same axis continuation (#370)
    !> as the Boozer field: a negative radial label s is the point
    !> (|s|, theta + pi), and the radial derivatives are those of the field
    !> as a function of the signed label. The symplectic integrators commit
    !> exactly this chart switch after a step through the axis.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: pi
    use simple, only: init_vmec
    use field, only: vmec_field_t, create_vmec_field
    use get_can_sub, only: get_canonical_coordinates_with_field
    use field_can_flux, only: evaluate_flux
    use field_can_base, only: field_can_t

    implicit none

    real(dp), parameter :: s = 0.02_dp, theta = 0.7_dp, phi = 0.3_dp
    real(dp), parameter :: ds = 1.0e-5_dp
    type(vmec_field_t) :: vmec_field
    type(field_can_t) :: fneg, fpos, fp, fm
    real(dp) :: field_period, fd

    call init_vmec('wout.nc', 5, 5, 3, field_period)
    call create_vmec_field(vmec_field)
    call get_canonical_coordinates_with_field(vmec_field)

    call evaluate_flux(fneg, -s, theta, phi, 2)
    call evaluate_flux(fpos, s, theta + pi, phi, 2)
    call check('Bmod', fneg%Bmod, fpos%Bmod)
    call check('Aph', fneg%Aph, fpos%Aph)
    call check('hth', fneg%hth, fpos%hth)
    call check('hph', fneg%hph, fpos%hph)
    call check('dBmod/ds', fneg%dBmod(1), -fpos%dBmod(1))
    call check('dBmod/dtheta', fneg%dBmod(2), fpos%dBmod(2))
    call check('dhph/ds', fneg%dhph(1), -fpos%dhph(1))
    call check('d2Bmod/ds2', fneg%d2Bmod(1), fpos%d2Bmod(1))
    call check('d2Bmod/dsdtheta', fneg%d2Bmod(2), -fpos%d2Bmod(2))

    ! The radial derivative at a negative label is the derivative along the
    ! signed label, which crosses into the opposite ray.
    call evaluate_flux(fp, -s + ds, theta, phi, 0)
    call evaluate_flux(fm, -s - ds, theta, phi, 0)
    fd = (fp%Bmod - fm%Bmod)/(2.0_dp*ds)
    if (abs(fd - fneg%dBmod(1)) > 1.0e-4_dp*abs(fneg%dBmod(1))) then
        print *, 'finite difference', fd, ' analytic', fneg%dBmod(1)
        error stop 'dBmod/ds at a negative label is not the signed derivative'
    end if
    print *, 'PASS'

contains

    subroutine check(name, actual, expected)
        character(*), intent(in) :: name
        real(dp), intent(in) :: actual, expected

        if (abs(actual - expected) > 1.0e-12_dp*max(1.0_dp, abs(expected))) then
            print *, name, ': at -s ', actual, ', at (s, theta+pi) ', expected
            error stop 'negative canonical-flux label is not the opposite ray'
        end if
    end subroutine check

end program test_canflux_axis_continuation
