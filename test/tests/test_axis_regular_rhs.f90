program test_axis_regular_rhs
    !> The RK recovery field in (X,Y,phi) must have a unique axis limit,
    !> preserve the established polar RHS outside its reconstruction disk,
    !> join without a handoff seam, and close after one field period.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: twopi
    use simple, only: tracer_t, init_params
    use simple_main, only: init_field
    use params, only: field_input, coord_input
    use velo_mod, only: isw_field_type
    use magfie_sub, only: BOOZER, init_magfie
    use alpha_lifetime_sub, only: velo_axis_regularized, velo_can, &
                                  axis_polar_to_cartesian, axis_cartesian_to_polar, &
                                  axis_native_radius_squared, axis_polar_velocity

    implicit none

    integer, parameter :: nray = 8
    real(dp), parameter :: sample_s = 1.0e-3_dp
    real(dp), parameter :: probe_rho = 1.0e-12_dp
    real(dp), parameter :: axis_limit_tol = 1.0e-8_dp
    type(tracer_t) :: norb
    real(dp) :: angle, spread, rho, s, theta
    real(dp) :: zaxis(5), vaxis(5), vprobe(5), zpolar(5), vpolar(5)
    real(dp) :: vexpected(5), vperiodic(5)
    integer :: iray
    logical :: failed

    isw_field_type = BOOZER
    field_input = 'wout.nc'
    coord_input = 'wout.nc'
    call init_field(norb, 'wout.nc', 5, 5, 5, 1)
    call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0e-12_dp)
    call init_magfie(BOOZER)

    failed = .false.
    call test_manufactured_axis_maps(failed)
    zaxis = [0.0_dp, 0.0_dp, 0.1_dp, 1.0_dp, -0.3_dp]
    call velo_axis_regularized(0.0_dp, zaxis, vaxis, sample_s)
    spread = 0.0_dp
    do iray = 0, nray - 1
        angle = twopi*real(iray, dp)/real(nray, dp)
        zaxis(1) = probe_rho*cos(angle)
        zaxis(2) = probe_rho*sin(angle)
        call velo_axis_regularized(0.0_dp, zaxis, vprobe, sample_s)
        spread = max(spread, maxval(abs(vprobe - vaxis)))
    end do
    if (spread > axis_limit_tol) then
        print *, 'FAIL: axis RHS depends on approach ray, spread = ', spread
        failed = .true.
    end if

    rho = 5.0e-2_dp
    theta = 0.7_dp
    s = rho*rho
    zaxis = [rho*cos(theta), rho*sin(theta), 0.1_dp, 1.0_dp, -0.3_dp]
    zpolar = [s, theta, zaxis(3), zaxis(4), zaxis(5)]
    call velo_can(0.0_dp, zpolar, vpolar)
    call polar_rhs_to_axis(zaxis, s, vpolar, vexpected)
    call velo_axis_regularized(0.0_dp, zaxis, vprobe, sample_s)
    if (maxval(abs(vprobe - vexpected)) > 1.0e-12_dp) then
        print *, 'FAIL: regular chart changed the outer RHS, error = ', &
            maxval(abs(vprobe - vexpected))
        failed = .true.
    end if

    rho = 0.999_dp*sqrt(sample_s)
    s = rho*rho
    zaxis(1:2) = [rho*cos(theta), rho*sin(theta)]
    zpolar(1:2) = [s, theta]
    call velo_can(0.0_dp, zpolar, vpolar)
    call polar_rhs_to_axis(zaxis, s, vpolar, vexpected)
    call velo_axis_regularized(0.0_dp, zaxis, vprobe, sample_s)
    if (maxval(abs(vprobe - vexpected)) > 1.0e-7_dp) then
        print *, 'FAIL: reconstructed axis RHS has a handoff seam, error = ', &
            maxval(abs(vprobe - vexpected))
        failed = .true.
    end if

    zaxis = [0.0_dp, 0.0_dp, 0.1_dp, 1.0_dp, -0.3_dp]
    call velo_axis_regularized(0.0_dp, zaxis, vaxis, sample_s)
    zaxis(3) = zaxis(3) + norb%fper
    call velo_axis_regularized(0.0_dp, zaxis, vperiodic, sample_s)
    if (maxval(abs(vperiodic - vaxis)) > 1.0e-10_dp) then
        print *, 'FAIL: axis RHS is not field-periodic, error = ', &
            maxval(abs(vperiodic - vaxis))
        failed = .true.
    end if

    if (failed) error stop 1
    print *, 'test_axis_regular_rhs PASSED'

contains

    subroutine test_manufactured_axis_maps(failed)
        logical, intent(inout) :: failed
        real(dp), parameter :: a = 0.37_dp, b = -0.22_dp, omega = 0.83_dp
        real(dp), parameter :: radii(5) = &
                               [1.0e-8_dp, 1.0e-5_dp, 1.0e-3_dp, 0.05_dp, 0.4_dp]
        real(dp), parameter :: angles(6) = &
                  [-0.9_dp, 0.0_dp, 0.7_dp, 0.25_dp*twopi, 0.5_dp*twopi, twopi + 0.3_dp]
        real(dp), parameter :: native_floors(4) = &
                               [0.0_dp, 1.0e-12_dp, 1.0e-4_dp, 0.2_dp]
        real(dp) :: polar(5), cart(5), recovered(5), expected(5), velocity(5)
        real(dp) :: expected_velocity(5), cart_velocity(5), recart(5)
        real(dp) :: radius, phase, radial_velocity, q, expected_radius2, x
        integer :: convention, ir, ia, j
        logical :: radial_is_rho

        do convention = 0, 1
            radial_is_rho = convention == 1
            do j = 1, size(native_floors)
                q = native_floors(j)
                expected_radius2 = q
                if (radial_is_rho) expected_radius2 = q*q
                if (abs(axis_native_radius_squared(q, radial_is_rho) - &
                        expected_radius2) > 1.0e-14_dp*max(expected_radius2, &
                                                           tiny(1.0_dp))) then
                    print *, 'FAIL: native radial floor has wrong squared units'
                    failed = .true.
                end if
            end do
            do ir = 1, size(radii)
                radius = radii(ir)
                do ia = 1, size(angles)
                    phase = angles(ia)
                    q = radius*radius
                    if (radial_is_rho) q = radius
                    polar = [q, phase, 0.19_dp, 0.8_dp, -0.4_dp]
                    expected = [radius*cos(phase), radius*sin(phase), polar(3:5)]
                    call axis_polar_to_cartesian(polar, cart, radial_is_rho)
                    if (maxval(abs(cart - expected)) > 2.0e-14_dp) then
                        print *, 'FAIL: manufactured polar-to-axis coordinate map'
                        failed = .true.
                    end if
                    call axis_cartesian_to_polar(expected, recovered, radial_is_rho)
                    if (abs(recovered(1) - q) > 2.0e-14_dp*q .or. &
                        abs(cos(recovered(2)) - cos(phase)) > 2.0e-14_dp .or. &
                        abs(sin(recovered(2)) - sin(phase)) > 2.0e-14_dp .or. &
                        any(recovered(3:5) /= polar(3:5))) then
                        print *, 'FAIL: manufactured axis-to-polar coordinate map'
                        failed = .true.
                    end if

                    ! Derive the polar field independently from Cartesian
                    ! translation plus rotation. Nonzero A and B expose the
                    ! missing radial chain factor that pure rotation hides.
                    radial_velocity = a*cos(phase) + b*sin(phase)
                    velocity(1) = 2.0_dp*radius*radial_velocity
                    if (radial_is_rho) velocity(1) = radial_velocity
                    velocity(2) = (b*cos(phase) - a*sin(phase))/radius + omega
                    velocity(3:5) = [0.17_dp, -0.03_dp, 0.04_dp]
                    expected_velocity = [a - omega*expected(2), &
                                         b + omega*expected(1), velocity(3:5)]
                    call axis_polar_velocity(expected, velocity, cart_velocity, &
                                             radial_is_rho)
                    if (maxval(abs(cart_velocity - expected_velocity)) > &
                        2.0e-12_dp) then
                        print *, 'FAIL: manufactured Cartesian RHS chain rule'
                        failed = .true.
                    end if
                end do
            end do

            ! A straight physical path crosses the axis without changing its
            ! Cartesian position or tail state. Theta itself is undefined at
            ! the origin, so test the round trip rather than choosing a ray.
            do j = -4, 4
                x = 0.025_dp*real(j, dp)
                expected = [x, 0.0_dp, 0.19_dp, 0.8_dp, -0.4_dp]
                call axis_cartesian_to_polar(expected, recovered, radial_is_rho)
                call axis_polar_to_cartesian(recovered, recart, radial_is_rho)
                if (maxval(abs(recart - expected)) > 2.0e-14_dp) then
                    print *, 'FAIL: manufactured straight axis crossing round trip'
                    failed = .true.
                end if
            end do
        end do
    end subroutine test_manufactured_axis_maps

    pure subroutine polar_rhs_to_axis(zaxis, s, vpolar, vaxis)
        real(dp), intent(in) :: zaxis(5), s, vpolar(5)
        real(dp), intent(out) :: vaxis(5)

        vaxis(1) = 0.5_dp*vpolar(1)*zaxis(1)/s - vpolar(2)*zaxis(2)
        vaxis(2) = 0.5_dp*vpolar(1)*zaxis(2)/s + vpolar(2)*zaxis(1)
        vaxis(3:5) = vpolar(3:5)
    end subroutine polar_rhs_to_axis

end program test_axis_regular_rhs
