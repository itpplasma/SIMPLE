program test_spectre_fo_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_spectre, only: spectre_field_t, create_spectre_field
    use magfie_sub, only: set_magfie_spectre_field, spectre_field, &
        TESLA_TO_GAUSS, M_TO_CM
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_enter, &
        spectre_fo_to_gc, &
        spectre_fo_canonical_pzeta, spectre_fo_advance_until_exit, &
        SPECTRE_FO_OK, SPECTRE_FO_LOSS
    implicit none

    real(dp), parameter :: RO0_BAR = 0.01_dp
    real(dp), parameter :: A_SCALE = TESLA_TO_GAUSS*M_TO_CM**2
    character(len=1024) :: filename
    type(spectre_field_t) :: field
    type(spectre_fo_state_t) :: state
    real(dp) :: y(5), y_back(5), acov(3), hcov(3), bmod, sqgb(3)
    real(dp) :: target, pzeta, entry_residual, energy0, speed0, x_exact(3)
    real(dp) :: pzeta_particle, pzeta_gc, mu_bar, dt_used, y_exit(5), u_eval(3)
    real(dp) :: e_cov(3, 3), g(3, 3), ginv(3, 3), sqrtg, bvec(3), radial(3)
    real(dp) :: gyro_period, pzeta_scale, probe_time
    real(dp) :: speed_residual
    integer :: ierr, owner
    logical :: exited

    if (command_argument_count() /= 1) error stop 'expected SPECTRE fixture'
    call get_command_argument(1, filename)
    call create_spectre_field(field, trim(filename), ierr)
    if (ierr /= 0) error stop 'failed to load SPECTRE fixture'
    call set_magfie_spectre_field(field)

    y = [1.0_dp, 0.7_dp, 0.3_dp, 0.5_dp, 0.7_dp]
    call spectre_field%evaluate([1.0_dp - 1.0e-12_dp, y(2), y(3)], &
        acov, hcov, bmod, sqgb)
    target = sqrt(2.0_dp)*y(4)*y(5)*hcov(3)*M_TO_CM + &
        acov(3)*A_SCALE/RO0_BAR
    call spectre_fo_enter(state, y, 1, 1, RO0_BAR, ierr)
    if (ierr /= SPECTRE_FO_OK) error stop 'GC to FO entry failed'
    energy0 = 0.5_dp*dot_product(state%v, state%v)
    if (abs(energy0 - y(4)**2) > 32.0_dp*epsilon(1.0_dp)*y(4)**2) &
        error stop 'entry energy mismatch'
    call spectre_fo_canonical_pzeta(state, RO0_BAR, pzeta, ierr)
    if (ierr /= SPECTRE_FO_OK) error stop 'P_zeta evaluation failed'
    entry_residual = abs(pzeta - target)/max(abs(target), 1.0_dp)
    if (entry_residual > 2.0e-10_dp) &
        error stop 'entry P_zeta mismatch'

    call spectre_field%coords%covariant_basis(state%u, e_cov)
    pzeta_scale = max(abs(dot_product(state%v, e_cov(:, 3)*M_TO_CM)), 1.0_dp)
    gyro_period = 8.0_dp*atan(1.0_dp)*RO0_BAR/(bmod*TESLA_TO_GAUSS)
    probe_time = 20.0_dp*gyro_period
    speed0 = norm2(state%v)
    call spectre_fo_advance_until_exit(state, probe_time, RO0_BAR, dt_used, &
        y_exit, owner, exited, ierr)
    if (ierr /= SPECTRE_FO_OK) error stop 'Boris step failed'
    if (exited) error stop 'premature full-orbit exit'
    if (abs(dt_used - probe_time) > epsilon(1.0_dp)*probe_time) &
        error stop 'full-orbit time accounting'
    speed_residual = abs(norm2(state%v) - speed0)/speed0
    if (speed_residual > 1.0e-12_dp) &
        error stop 'Boris speed changed'

    call spectre_fo_canonical_pzeta(state, RO0_BAR, pzeta_particle, ierr)
    if (ierr /= SPECTRE_FO_OK) error stop 'particle P_zeta evaluation failed'
    if (abs(pzeta_particle - target) > 1.0e-2_dp*pzeta_scale) &
        error stop 'Boris P_zeta drift'
    call spectre_fo_to_gc(state, RO0_BAR, y_back, owner, ierr, mu_bar)
    if (ierr /= SPECTRE_FO_OK) error stop 'FO to GC reconstruction failed'
    if (abs(y_back(4)**2 - energy0) > 1.0e-12_dp*energy0) &
        error stop 'reconstruction energy mismatch'
    u_eval = y_back(1:3)
    if (owner < spectre_field%data%Mvol) then
        if (u_eval(1) == real(owner, dp)) u_eval(1) = u_eval(1) - 1.0e-12_dp
    end if
    call spectre_field%evaluate(u_eval, acov, hcov, bmod, sqgb)
    pzeta_gc = sqrt(2.0_dp)*y_back(4)*y_back(5)*hcov(3)*M_TO_CM + &
        acov(3)*A_SCALE/RO0_BAR
    if (abs(pzeta_gc - pzeta_particle) > &
        2.0e-10_dp*max(abs(pzeta_particle), 1.0_dp)) &
        error stop 'reconstruction P_zeta mismatch'
    if (abs(mu_bar - y_back(4)**2*(1.0_dp - y_back(5)**2)/ &
        (bmod*TESLA_TO_GAUSS)) > 64.0_dp*epsilon(1.0_dp)*max(mu_bar, 1.0_dp)) &
        error stop 'reconstruction mu mismatch'

    call spectre_field%coords%evaluate_cart([1.0_dp, y(2), y(3)], x_exact)
    state%x = x_exact*M_TO_CM
    state%u = [1.0_dp, y(2), y(3)]
    state%owner = 1
    call spectre_fo_canonical_pzeta(state, RO0_BAR, pzeta, ierr)
    if (ierr /= SPECTRE_FO_OK) error stop 'exact-interface locate failed'
    if (state%owner /= 1) error stop 'exact interface changed owner'

    state = spectre_fo_state_t()
    state%u = [1.5_dp, y(2), y(3)]
    state%owner = 2
    state%iface = 1
    state%active = .true.
    call spectre_field%coords%evaluate_cart(state%u, state%x)
    state%x = state%x*M_TO_CM
    call spectre_field%evaluate(state%u, acov, hcov, bmod, sqgb)
    call spectre_field%coords%covariant_basis(state%u, e_cov)
    call spectre_field%coords%metric_tensor(state%u, g, ginv, sqrtg)
    bvec = matmul(e_cov, sqgb/sqrtg)
    state%v = bvec/norm2(bvec)
    call spectre_fo_advance_until_exit(state, 1.0e-8_dp, RO0_BAR, dt_used, &
        y_exit, owner, exited, ierr)
    if (ierr /= SPECTRE_FO_OK .or. .not. exited .or. owner /= 2) &
        error stop 'full-orbit return to GC failed'
    if (dt_used <= 0.0_dp .or. dt_used > 1.0e-8_dp) &
        error stop 'full-orbit exit time accounting'

    state%u = [real(spectre_field%data%Mvol, dp), y(2), y(3)]
    state%owner = spectre_field%data%Mvol
    state%iface = 1
    state%active = .true.
    call spectre_field%coords%evaluate_cart(state%u, state%x)
    state%x = state%x*M_TO_CM
    call spectre_field%coords%covariant_basis(state%u, e_cov)
    radial = e_cov(:, 1)/norm2(e_cov(:, 1))
    state%v = radial
    call spectre_fo_advance_until_exit(state, 1.0e-4_dp, RO0_BAR, dt_used, &
        y_exit, owner, exited, ierr)
    if (ierr /= SPECTRE_FO_LOSS .or. &
        state%u(1) /= real(state%owner, dp)) &
        error stop 'outer particle boundary loss failed'
    call spectre_field%coords%evaluate_cart(state%u, state%x)
    state%x = state%x*M_TO_CM
    state%v = -radial
    state%active = .true.
    call spectre_fo_advance_until_exit(state, 1.0e-8_dp, RO0_BAR, dt_used, &
        y_exit, owner, exited, ierr)
    if (ierr /= SPECTRE_FO_OK) error stop 'inward outer-boundary state lost'

    print '(A,ES12.4)', 'spectre_fo: entry relative P_zeta residual = ', &
        entry_residual
    print '(A,ES12.4)', 'spectre_fo: relative speed residual = ', &
        speed_residual
end program test_spectre_fo_hybrid
