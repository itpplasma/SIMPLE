program test_rk_angle_seam_crossing
    !> Regression test for adaptive RK guiding-centre steps in VMEC
    !> coordinates whose poloidal or toroidal angle passes through zero.
    !>
    !> The VMEC splines are continuous at their knots, including 0 = 2*pi,
    !> only to about 1e-8 relative in the right-hand side near the boundary
    !> of the public QA equilibrium. With a purely relative error scale, the
    !> scale of an angle crossing zero shrinks like the step itself, so the
    !> error estimate of a step over the jump no longer decreases with the
    !> step and the step size underflows (unresolved markers in vmec_rk runs).
    !> The step must succeed and agree with the same step taken one period
    !> further on.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use magfie_sub, only: init_magfie, VMEC
    use alpha_lifetime_sub, only: orbit_timestep_axis, velo_can
    use parmot_mod, only: ro0, rmu

    implicit none

    real(dp), parameter :: dtau = 25.0_dp, relerr = 1.0e-13_dp
    integer, parameter :: n_start = 4
    real(dp) :: fper, twopi, s0(n_start), angle0(n_start), lambda0(n_start)
    real(dp) :: z0(5), vz(5), z_seam(5), z_shift(5), deviation
    integer :: i, iangle, ierr_seam, ierr_shift
    logical :: failed

    call init_vmec('wout.nc', 5, 5, 5, fper)
    ! Larmor-radius factor of a 3.5 MeV alpha particle [G cm], as in params.
    ro0 = 2.7e5_dp
    rmu = 1.0e8_dp
    call init_magfie(VMEC)

    twopi = 2.0_dp*acos(-1.0_dp)
    s0 = [0.95_dp, 0.98_dp, 0.99_dp, 0.996_dp]
    angle0 = [0.4_dp, 1.3_dp, 2.2_dp, 0.9_dp]
    lambda0 = [-0.36_dp, -0.33_dp, 0.3_dp, -0.36_dp]

    failed = .false.
    do iangle = 2, 3
        do i = 1, n_start
            ! Start half a step before the seam of angle iangle.
            z0 = [s0(i), angle0(i), angle0(i), 1.0_dp, lambda0(i)]
            z0(iangle) = 0.0_dp
            call velo_can(0.0_dp, z0, vz)
            z0(iangle) = -0.5_dp*dtau*vz(iangle)

            z_seam = z0
            call orbit_timestep_axis(z_seam, dtau, dtau, relerr, ierr_seam)
            ! The same start one period on meets the seam at |angle| = 2*pi.
            z_shift = z0
            z_shift(iangle) = z_shift(iangle) + sign(twopi, -vz(iangle))
            call orbit_timestep_axis(z_shift, dtau, dtau, relerr, ierr_shift)
            z_shift(iangle) = z_shift(iangle) - sign(twopi, -vz(iangle))

            deviation = maxval(abs(z_seam - z_shift))
            print '(a,i2,a,i2,a,2i3,a,es12.4)', 'angle ', iangle, ' start ', i, &
                ' ierr ', ierr_seam, ierr_shift, ' deviation ', deviation
            if (ierr_seam /= 0 .or. ierr_shift /= 0) then
                failed = .true.
            else if (deviation > 1.0e-9_dp) then
                failed = .true.
            else if (z_seam(iangle)*z0(iangle) >= 0.0_dp) then
                print *, 'step did not cross the seam'
                failed = .true.
            end if
        end do
    end do

    if (failed) error stop 'RK step across an angular seam failed'
    print *, 'PASSED: RK steps across theta = 0 and phi = 0'

end program test_rk_angle_seam_crossing
