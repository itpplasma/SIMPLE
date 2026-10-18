program test_rk_pitch_zero_start
    !> Regression test for adaptive RK guiding-centre steps started at zero
    !> pitch where the parallel derivative of B vanishes.
    !>
    !> There d(lambda)/dtau = 0 at lambda = 0, so lambda stays at roundoff size.
    !> A purely relative error scale then forces the step size towards zero and
    !> the step limit is exhausted (native exit code 104 in RK runs). The step
    !> must instead succeed and agree with the neighbouring starts lambda=+-delta.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use magfie_sub, only: init_magfie, VMEC
    use alpha_lifetime_sub, only: orbit_timestep_axis, velo_can
    use parmot_mod, only: ro0, rmu

    implicit none

    real(dp), parameter :: s0 = 0.25_dp, phi0 = 0.1_dp
    real(dp), parameter :: dtau = 100.0_dp, relerr = 1.0e-10_dp
    real(dp), parameter :: delta = 1.0e-10_dp
    integer, parameter :: step_limit = 20000, n_scan = 64

    real(dp) :: fper, theta_a, theta_b, g_a, g_b, theta_root
    real(dp) :: z0(5), z_zero(5), z_plus(5), z_minus(5), deviation
    integer :: i, ierr_zero, ierr_plus, ierr_minus

    call init_vmec('wout.nc', 5, 5, 5, fper)
    ! Larmor-radius factor of a 3.5 MeV alpha particle [G cm], as in params.
    ro0 = 2.7e5_dp
    rmu = 1.0e8_dp
    call init_magfie(VMEC)

    ! Bracket a zero of d(lambda)/dtau at lambda=0 in theta at fixed s and phi.
    theta_a = 0.0_dp
    g_a = pitch_rate(theta_a)
    theta_b = theta_a
    g_b = g_a
    do i = 1, n_scan
        theta_b = 2.0_dp*acos(-1.0_dp)*real(i, dp)/real(n_scan, dp)
        g_b = pitch_rate(theta_b)
        if (g_a*g_b <= 0.0_dp) exit
        theta_a = theta_b
        g_a = g_b
    end do
    if (g_a*g_b > 0.0_dp) error stop 'no zero of the mirror force found'

    do i = 1, 200
        theta_root = 0.5_dp*(theta_a + theta_b)
        if (theta_root <= theta_a .or. theta_root >= theta_b) exit
        if (g_a*pitch_rate(theta_root) <= 0.0_dp) then
            theta_b = theta_root
        else
            theta_a = theta_root
            g_a = pitch_rate(theta_a)
        end if
    end do
    print '(a,es24.16,a,es12.4)', 'theta at mirror-force zero: ', theta_root, &
        ' dlambda/dtau: ', pitch_rate(theta_root)

    z0 = [s0, theta_root, phi0, 1.0_dp, 0.0_dp]
    z_zero = z0
    call orbit_timestep_axis(z_zero, dtau, dtau, relerr, ierr_zero, step_limit)
    z_plus = z0
    z_plus(5) = delta
    call orbit_timestep_axis(z_plus, dtau, dtau, relerr, ierr_plus, step_limit)
    z_minus = z0
    z_minus(5) = -delta
    call orbit_timestep_axis(z_minus, dtau, dtau, relerr, ierr_minus, step_limit)

    print '(a,3i3)', 'ierr (lambda=0, +delta, -delta): ', ierr_zero, ierr_plus, &
        ierr_minus
    print '(a,5es24.16)', 'z(lambda=0):  ', z_zero
    print '(a,5es24.16)', 'z(+delta):    ', z_plus
    if (ierr_zero /= 0) error stop 'RK step from zero pitch failed'
    if (ierr_plus /= 0 .or. ierr_minus /= 0) &
        error stop 'RK step from neighbouring pitch failed'

    ! The solution is Lipschitz in lambda_0; the zero start must lie between
    ! (within tolerance of) its neighbours rather than being frozen in place.
    deviation = maxval(abs(z_zero(1:4) - 0.5_dp*(z_plus(1:4) + z_minus(1:4))))
    deviation = max(deviation, abs(z_zero(5) - 0.5_dp*(z_plus(5) + z_minus(5))))
    print '(a,es12.4)', 'deviation from neighbour mean: ', deviation
    if (deviation > 1.0e-7_dp) error stop 'zero-pitch step inconsistent with neighbours'
    if (all(z_zero(1:3) == z0(1:3))) error stop 'zero-pitch state did not advance'

    print *, 'PASSED: RK step from zero pitch at a mirror-force zero'

contains

    real(dp) function pitch_rate(theta)
        real(dp), intent(in) :: theta
        real(dp) :: z(5), vz(5)

        z = [s0, theta, phi0, 1.0_dp, 0.0_dp]
        call velo_can(0.0_dp, z, vz)
        pitch_rate = vz(5)
    end function pitch_rate

end program test_rk_pitch_zero_start
