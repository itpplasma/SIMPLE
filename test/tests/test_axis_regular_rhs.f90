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
    use alpha_lifetime_sub, only: velo_axis_regularized, velo_can

    implicit none

    integer, parameter :: nray = 8
    real(dp), parameter :: sample_s = 1.0e-6_dp
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

    rho = 2.0e-2_dp
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

    pure subroutine polar_rhs_to_axis(zaxis, s, vpolar, vaxis)
        real(dp), intent(in) :: zaxis(5), s, vpolar(5)
        real(dp), intent(out) :: vaxis(5)

        vaxis(1) = 0.5_dp*vpolar(1)*zaxis(1)/s - vpolar(2)*zaxis(2)
        vaxis(2) = 0.5_dp*vpolar(1)*zaxis(2)/s + vpolar(2)*zaxis(1)
        vaxis(3:5) = vpolar(3:5)
    end subroutine polar_rhs_to_axis

end program test_axis_regular_rhs
