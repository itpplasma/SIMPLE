module linear_radial_field_backend
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use field_can_base, only: field_can_t
  use orbit_symplectic_base, only: symplectic_integrator_t, &
    SYMPLECTIC_STEP_BOUNDARY_LIMITED, SYMPLECTIC_STEP_OK, &
    SYMPLECTIC_STEP_OUTSIDE_DOMAIN, SYMPLECTIC_STEP_MAXITER, &
    SYMPLECTIC_STEP_BOUNDARY

  implicit none

  integer :: retry_calls = 0

contains

  subroutine evaluate_linear_radial(f, r, theta, phi, mode)
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: r, theta, phi
    integer, intent(in) :: mode

    f%Ath = r
    f%Aph = 0.0_dp
    f%hth = 0.0_dp
    f%hph = 1.0_dp
    f%Bmod = 2.0_dp - theta
    f%dAth = [1.0_dp, 0.0_dp, 0.0_dp]
    f%dAph = 0.0_dp
    f%dhth = 0.0_dp
    f%dhph = 0.0_dp
    f%dBmod = [0.0_dp, -1.0_dp, 0.0_dp]
    f%d2Ath = 0.0_dp
    f%d2Aph = 0.0_dp
    f%d2hth = 0.0_dp
    f%d2hph = 0.0_dp
    f%d2Bmod = 0.0_dp
  end subroutine evaluate_linear_radial

  subroutine basin_limited_step(si, f, step_status)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: step_status

    if (si%dt >= 0.1_dp) then
      step_status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
    else if (si%dt >= 0.0625_dp) then
      step_status = SYMPLECTIC_STEP_BOUNDARY_LIMITED
    else
      si%z(1) = 0.5_dp + (0.5_dp - 1.0e-12_dp)* &
        (1.0_dp - ((0.0625_dp - si%dt)/0.0625_dp)**2)
      step_status = SYMPLECTIC_STEP_OK
    end if
  end subroutine basin_limited_step

  subroutine nonlinear_boundary_step(si, f, step_status)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: step_status
    real(dp) :: radius

    if (si%dt > 0.4_dp) then
      step_status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
    else
      radius = 1.0_dp - 1.0e-7_dp - (0.4_dp - si%dt)**2
      si%z(1) = radius
      step_status = SYMPLECTIC_STEP_OK
    end if
  end subroutine nonlinear_boundary_step

  subroutine retryable_step(si, f, step_status)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: step_status

    retry_calls = retry_calls + 1
    if (si%dt > 0.5_dp) then
      si%z = 99.0_dp
      f%H = 99.0_dp
      step_status = SYMPLECTIC_STEP_MAXITER
      return
    end if
    si%z(1) = si%z(1) + si%dt
    f%H = f%H + si%dt
    step_status = SYMPLECTIC_STEP_OK
  end subroutine retryable_step

  subroutine failed_boundary_step(si, f, step_status)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: step_status

    retry_calls = retry_calls + 1
    si%z = 99.0_dp
    f%H = 99.0_dp
    step_status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
  end subroutine failed_boundary_step

  subroutine second_half_fails_step(si, f, step_status)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: step_status

    retry_calls = retry_calls + 1
    if (si%dt > 0.5_dp .or. si%z(1) > 0.5_dp) then
      step_status = SYMPLECTIC_STEP_MAXITER
      return
    end if
    si%z(1) = si%z(1) + si%dt
    f%H = f%H + si%dt
    step_status = SYMPLECTIC_STEP_OK
  end subroutine second_half_fails_step

  subroutine retryable_boundary_step(si, f, step_status)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: step_status

    retry_calls = retry_calls + 1
    if (si%dt > 0.5_dp) then
      step_status = SYMPLECTIC_STEP_MAXITER
      return
    end if
    if (retry_calls == 2) then
      si%z(1) = si%z(1) + si%dt
      f%H = f%H + si%dt
      step_status = SYMPLECTIC_STEP_OK
      return
    end if
    si%z(1) = si%z(1) + 0.25_dp*si%dt
    f%H = f%H + 0.25_dp*si%dt
    si%last_step_fraction = 0.25_dp
    si%last_event_fraction_width = 0.02_dp
    step_status = SYMPLECTIC_STEP_BOUNDARY
  end subroutine retryable_boundary_step
end module linear_radial_field_backend

program test_newton_solver_status
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use field_can_mod, only: field_can_t, eval_field => evaluate
  use linear_radial_field_backend, only: basin_limited_step, &
    evaluate_linear_radial, failed_boundary_step, nonlinear_boundary_step, &
    retryable_boundary_step, retryable_step, retry_calls, second_half_fails_step
  use orbit_symplectic, only: guard_lobatto_stage_radii, boundary_event_converged, &
    advance_symplectic_with_boundary, advance_symplectic_with_retry, &
    newton_midpoint, orbit_sympl_init, &
    orbit_timestep_sympl, matrix3_near_singular, solve_newton_system, &
    get_boundary_event_tolerances, accept_warning_maxiter
  use orbit_symplectic_base, only: symplectic_integrator_t, &
    EXPL_IMPL_EULER, IMPL_EXPL_EULER, MIDPOINT, GAUSS1, GAUSS2, GAUSS3, &
    GAUSS4, LOBATTO3, &
    SYMPLECTIC_STEP_BOUNDARY, &
    SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_OUTSIDE_DOMAIN, &
    SYMPLECTIC_STEP_MAXITER, SYMPLECTIC_STEP_LINEAR_SOLVE, &
    SYMPLECTIC_STEP_EVENT_NOT_CONVERGED, SYMPLECTIC_STEP_BOUNDARY_LIMITED, &
    boundary_event_fraction_tolerance, boundary_event_radial_tolerance, &
    symplectic_newton_warning_mode

  implicit none

  type(symplectic_integrator_t) :: integrator
  type(field_can_t) :: field
  real(dp) :: x(5), xlast(5), matrix(2, 2), residual(2), matrix3(3, 3)
  real(dp) :: lobatto_state(10), expected_lobatto_state(10)
  integer :: i, status

  field%Aph = 0.0_dp
  field%ro0 = 1.0_dp
  x = [0.5_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.5_dp]
  symplectic_newton_warning_mode = .false.
  call newton_midpoint(integrator, field, x, 1.0e-15_dp, 1.0e-12_dp, &
    0, xlast, status)
  if (status /= SYMPLECTIC_STEP_MAXITER) then
    error stop 'zero-iteration midpoint solve did not report max iterations'
  end if
  if (any(xlast /= x)) error stop 'zero-iteration solve changed the iterate'

  x(1) = ieee_value(0.0_dp, ieee_quiet_nan)
  symplectic_newton_warning_mode = .true.
  call newton_midpoint(integrator, field, x, 1.0e-15_dp, 1.0e-12_dp, &
    0, xlast, status)
  if (status /= SYMPLECTIC_STEP_MAXITER) then
    error stop 'warning mode accepted a non-finite Newton iterate'
  end if
  symplectic_newton_warning_mode = .false.

  x(1) = 1.1_dp
  call newton_midpoint(integrator, field, x, 1.0e-15_dp, 1.0e-12_dp, &
    1, xlast, status)
  if (status /= SYMPLECTIC_STEP_OUTSIDE_DOMAIN) then
    error stop 'exterior midpoint iterate did not report its domain failure'
  end if

  matrix = 0.0_dp
  residual = 1.0_dp
  call solve_newton_system(matrix, residual, status)
  if (status /= SYMPLECTIC_STEP_LINEAR_SOLVE) then
    error stop 'singular Newton system did not report its LAPACK failure'
  end if

  matrix = 0.0_dp
  matrix(1, 1) = 1.0_dp
  matrix(2, 2) = 1.0_dp
  residual = [2.0_dp, -3.0_dp]
  call solve_newton_system(matrix, residual, status)
  if (status /= SYMPLECTIC_STEP_OK) then
    error stop 'nonsingular Newton system reported a failure'
  end if
  if (any(residual /= [2.0_dp, -3.0_dp])) then
    error stop 'Newton system returned the wrong correction'
  end if

  matrix3 = 0.0_dp
  matrix3(1, 1) = 1.0e-8_dp
  matrix3(2, 2) = 1.0e-8_dp
  matrix3(3, 3) = 1.0e-8_dp
  if (matrix3_near_singular(matrix3, 1.0e-24_dp)) then
    error stop 'scaled well-conditioned Euler system was classified singular'
  end if

  lobatto_state = [(-0.1_dp*real(i, dp), i = 1, 10)]
  lobatto_state(1) = -0.2_dp
  lobatto_state(3) = -0.4_dp
  lobatto_state(7) = 0.7_dp
  expected_lobatto_state = lobatto_state
  expected_lobatto_state(1) = 0.01_dp
  expected_lobatto_state(3) = 0.01_dp
  call guard_lobatto_stage_radii(lobatto_state, 3, status)
  if (status /= SYMPLECTIC_STEP_OK) then
    error stop 'valid Lobatto stage radii reported a domain failure'
  end if
  if (any(lobatto_state /= expected_lobatto_state)) then
    error stop 'Lobatto radius guard changed a nonradial stage component'
  end if

  lobatto_state(7) = 1.1_dp
  call guard_lobatto_stage_radii(lobatto_state, 3, status)
  if (status /= SYMPLECTIC_STEP_OUTSIDE_DOMAIN) then
    error stop 'exterior Lobatto stage did not report its domain failure'
  end if

  if (.not. boundary_event_converged(1.0e-11_dp, 2.0e-9_dp, &
      1.0e-10_dp, 1.0e-8_dp)) then
    error stop 'converged boundary bracket was rejected'
  end if
  if (boundary_event_converged(1.0e-11_dp, 2.0e-5_dp, &
      1.0e-10_dp, 1.0e-8_dp)) then
    error stop 'large boundary radial residual was accepted'
  end if
  if (boundary_event_converged(1.0e-4_dp, 2.0e-9_dp, &
      1.0e-10_dp, 1.0e-8_dp)) then
    error stop 'wide boundary time bracket was accepted'
  end if

  call test_lcfs_location
  call test_configured_event_tolerances
  call test_solver_basin_is_not_boundary
  call test_newton_warning_mode
  call test_step_retry

contains

  subroutine test_newton_warning_mode
    integer, parameter :: modes(8) = [EXPL_IMPL_EULER, IMPL_EXPL_EULER, &
      MIDPOINT, GAUSS1, GAUSS2, GAUSS3, GAUSS4, LOBATTO3]
    type(symplectic_integrator_t) :: strict_integrator
    type(field_can_t) :: strict_field
    real(dp) :: initial_state(4), accepted(2), previous(2), tolref(2)
    integer :: mode_index, step_status

    eval_field => evaluate_linear_radial
    initial_state = [0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    do mode_index = 1, size(modes)
      strict_field%ro0 = 1.0_dp
      strict_field%mu = 1.0_dp
      call orbit_sympl_init(strict_integrator, strict_field, initial_state, &
        0.01_dp, 1, 0.0_dp, modes(mode_index))
      strict_integrator%atol = 0.0_dp
      symplectic_newton_warning_mode = .false.
      call orbit_timestep_sympl(strict_integrator, strict_field, step_status)
      if (step_status /= SYMPLECTIC_STEP_MAXITER) then
        print *, 'strict mode, status:', modes(mode_index), step_status
        error stop 'strict timestep did not report max iterations'
      end if
      if (any(strict_integrator%z /= initial_state)) then
        error stop 'strict max-iteration failure changed the accepted state'
      end if

    end do

    previous = [1.0_dp, 2.0_dp]
    tolref = [1.0_dp, 2.0_dp]
    accepted = previous + [5.0e-12_dp, 1.0e-11_dp]
    symplectic_newton_warning_mode = .true.
    if (.not. accept_warning_maxiter(accepted, previous, tolref, 1.0e-12_dp)) then
      error stop 'warning mode rejected a bounded Newton correction'
    end if
    accepted(1) = huge(1.0_dp)
    if (accept_warning_maxiter(accepted, previous, tolref, 1.0e-12_dp)) then
      error stop 'warning mode accepted an unbounded Newton correction'
    end if
    accepted = previous
    accepted(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    if (accept_warning_maxiter(accepted, previous, tolref, 1.0e-12_dp)) then
      error stop 'warning mode accepted a non-finite Newton correction'
    end if
    accepted = previous
    symplectic_newton_warning_mode = .false.
    if (accept_warning_maxiter(accepted, previous, tolref, 1.0e-12_dp)) then
      error stop 'strict mode accepted a Newton max-iteration state'
    end if
  end subroutine test_newton_warning_mode

  subroutine test_step_retry
    type(symplectic_integrator_t) :: retry_integrator
    type(field_can_t) :: retry_field
    real(dp) :: initial_state(4), accepted_fraction
    integer :: step_status

    initial_state = [0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    retry_integrator%z = initial_state
    retry_integrator%dt = 1.0_dp
    retry_field%H = 0.0_dp
    retry_calls = 0
    symplectic_newton_warning_mode = .true.
    call advance_symplectic_with_retry(retry_integrator, retry_field, &
      retryable_step, step_status)
    if (step_status /= SYMPLECTIC_STEP_OK) then
      error stop 'warning mode did not recover a failed full step'
    end if
    if (retry_calls /= 3) error stop 'warning mode used the wrong retry sequence'
    if (abs(retry_integrator%z(1) - 1.5_dp) > 1.0e-14_dp) then
      error stop 'recovered step did not advance the full interval'
    end if
    if (abs(retry_field%H - 1.0_dp) > 1.0e-14_dp) then
      error stop 'recovered field state did not advance the full interval'
    end if
    if (retry_integrator%dt /= 1.0_dp) then
      error stop 'recovered step did not restore the configured timestep'
    end if

    retry_integrator%z = initial_state
    retry_integrator%dt = 1.0_dp
    retry_field%H = 0.0_dp
    retry_calls = 0
    call advance_symplectic_with_retry(retry_integrator, retry_field, &
      second_half_fails_step, step_status, accepted_fraction)
    if (step_status /= SYMPLECTIC_STEP_OK) then
      error stop 'warning mode discarded a recoverable first half'
    end if
    if (abs(retry_integrator%z(1) - 1.0_dp) > 1.0e-14_dp .or. &
        abs(retry_field%H - 0.5_dp) > 1.0e-14_dp) then
      error stop 'failed second half rolled back accepted progress'
    end if
    if (retry_integrator%dt /= 1.0_dp) then
      error stop 'partial retry did not restore the configured timestep'
    end if
    if (abs(accepted_fraction - 0.5_dp) > 1.0e-14_dp) then
      error stop 'partial retry reported the wrong accepted duration'
    end if

    retry_integrator%z = initial_state
    retry_integrator%dt = 1.0_dp
    retry_field%H = 0.0_dp
    retry_calls = 0
    symplectic_newton_warning_mode = .false.
    call advance_symplectic_with_retry(retry_integrator, retry_field, &
      retryable_step, step_status)
    if (step_status /= SYMPLECTIC_STEP_MAXITER .or. retry_calls /= 1) then
      error stop 'strict mode retried a failed step'
    end if
    if (any(retry_integrator%z /= initial_state) .or. retry_field%H /= 0.0_dp) then
      error stop 'failed strict step changed the accepted state'
    end if

    retry_calls = 0
    symplectic_newton_warning_mode = .true.
    call advance_symplectic_with_retry(retry_integrator, retry_field, &
      failed_boundary_step, step_status)
    if (step_status /= SYMPLECTIC_STEP_OUTSIDE_DOMAIN .or. retry_calls /= 1) then
      error stop 'warning mode retried a physical boundary status'
    end if
    if (any(retry_integrator%z /= initial_state) .or. retry_field%H /= 0.0_dp) then
      error stop 'failed boundary step changed the accepted state'
    end if

    retry_calls = 0
    call advance_symplectic_with_retry(retry_integrator, retry_field, &
      retryable_boundary_step, step_status)
    if (step_status /= SYMPLECTIC_STEP_BOUNDARY .or. retry_calls /= 3) then
      error stop 'retry path lost a converged boundary event'
    end if
    if (abs(retry_integrator%z(1) - 1.125_dp) > 1.0e-14_dp .or. &
        abs(retry_field%H - 0.625_dp) > 1.0e-14_dp) then
      error stop 'retry path lost the boundary event state'
    end if
    if (abs(retry_integrator%last_step_fraction - 0.625_dp) > 1.0e-14_dp .or. &
        abs(retry_integrator%last_event_fraction_width - 0.01_dp) > &
        1.0e-14_dp) then
      error stop 'retry path reported the wrong boundary event time'
    end if
  end subroutine test_step_retry

  subroutine test_configured_event_tolerances
    type(symplectic_integrator_t) :: loose_integrator, tight_integrator
    type(field_can_t) :: event_field
    real(dp) :: configured_fraction_tolerance, configured_radial_tolerance
    integer :: event_status

    eval_field => evaluate_linear_radial
    event_field%ro0 = 1.0_dp
    event_field%mu = 1.0_dp

    loose_integrator%z = [0.9_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    loose_integrator%dt = 1.0_dp
    loose_integrator%ntau = 1
    loose_integrator%rtol = 1.0e-12_dp
    boundary_event_fraction_tolerance = 1.0e-4_dp
    boundary_event_radial_tolerance = 1.0e-4_dp
    call get_boundary_event_tolerances(1.0e-12_dp, &
      configured_fraction_tolerance, configured_radial_tolerance)
    if (configured_fraction_tolerance /= 1.0e-4_dp .or. &
        configured_radial_tolerance /= 1.0e-4_dp) then
      error stop 'configured event tolerances were not selected'
    end if
    call advance_symplectic_with_boundary(loose_integrator, event_field, &
      nonlinear_boundary_step, event_status)
    if (event_status /= SYMPLECTIC_STEP_BOUNDARY) then
      error stop 'loose configured event tolerances rejected the crossing'
    end if

    tight_integrator%z = [0.9_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    tight_integrator%dt = 1.0_dp
    tight_integrator%ntau = 1
    tight_integrator%rtol = 1.0e-12_dp
    boundary_event_fraction_tolerance = 1.0e-8_dp
    boundary_event_radial_tolerance = 1.0e-8_dp
    call advance_symplectic_with_boundary(tight_integrator, event_field, &
      nonlinear_boundary_step, event_status)
    boundary_event_fraction_tolerance = -1.0_dp
    boundary_event_radial_tolerance = -1.0_dp
    if (event_status /= SYMPLECTIC_STEP_EVENT_NOT_CONVERGED) then
      error stop 'tight radial tolerance accepted an unresolved crossing'
    end if
    if (any(tight_integrator%z /= [0.9_dp, 0.0_dp, 0.0_dp, 0.0_dp])) then
      error stop 'unresolved configured event changed the accepted state'
    end if
  end subroutine test_configured_event_tolerances

  subroutine test_solver_basin_is_not_boundary
    type(symplectic_integrator_t) :: basin_integrator
    type(field_can_t) :: basin_field
    real(dp) :: initial_state(4)
    integer :: basin_status

    initial_state = [0.5_dp, 0.1_dp, 0.2_dp, 0.3_dp]
    basin_integrator%z = initial_state
    basin_integrator%dt = 0.1_dp
    basin_integrator%ntau = 1
    basin_integrator%rtol = 1.0e-12_dp
    call advance_symplectic_with_boundary(basin_integrator, basin_field, &
      basin_limited_step, basin_status)
    if (basin_status /= SYMPLECTIC_STEP_EVENT_NOT_CONVERGED) then
      error stop 'solver-basin limit was promoted to a physical boundary event'
    end if
    if (any(basin_integrator%z /= initial_state)) then
      error stop 'solver-basin failure changed the accepted state'
    end if
  end subroutine test_solver_basin_is_not_boundary

  subroutine test_lcfs_location
    real(dp), parameter :: timesteps(3) = [0.02_dp, 0.01_dp, 0.005_dp]
    integer, parameter :: modes(4) = [EXPL_IMPL_EULER, IMPL_EXPL_EULER, &
      MIDPOINT, LOBATTO3]
    real(dp), parameter :: exact_event_time = 0.1_dp
    real(dp) :: event_time, current_error, previous_error
    integer :: i, j

    do j = 1, size(modes)
      previous_error = huge(1.0_dp)
      do i = 1, size(timesteps)
        call locate_linear_boundary(modes(j), timesteps(i), event_time)
        current_error = abs(event_time - exact_event_time)
        if (current_error > 2.0_dp*timesteps(i)) then
          print *, 'mode, timestep, event time:', modes(j), timesteps(i), &
            event_time
          error stop 'boundary time disagrees with analytic crossing'
        end if
        if (current_error > previous_error + 256.0_dp*epsilon(1.0_dp)) then
          print *, 'mode, timestep, errors:', modes(j), timesteps(i), &
            previous_error, current_error
          error stop 'boundary time failed to converge under timestep refinement'
        end if
        previous_error = current_error
      end do
    end do
  end subroutine test_lcfs_location

  subroutine locate_linear_boundary(mode, timestep, event_time)
    integer, intent(in) :: mode
    real(dp), intent(in) :: timestep
    real(dp), intent(out) :: event_time
    type(symplectic_integrator_t) :: event_integrator
    type(field_can_t) :: event_field
    real(dp) :: z0(4)
    integer :: event_status, step

    eval_field => evaluate_linear_radial
    event_field%ro0 = 1.0_dp
    event_field%mu = 1.0_dp
    z0 = [0.9_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call orbit_sympl_init(event_integrator, event_field, z0, &
      timestep, 1, 1.0e-12_dp, mode)

    do step = 1, 1000
      call orbit_timestep_sympl(event_integrator, event_field, event_status)
      if (event_status == SYMPLECTIC_STEP_BOUNDARY) exit
      if (event_status /= SYMPLECTIC_STEP_OK) then
        print *, 'mode, timestep, status, residual, width:', mode, timestep, &
          event_status, event_integrator%last_event_radial_residual, &
          event_integrator%last_event_fraction_width
        error stop 'analytic boundary trace ended with a numerical failure'
      end if
    end do
    if (event_status /= SYMPLECTIC_STEP_BOUNDARY) then
      error stop 'analytic boundary trace did not reach the boundary'
    end if
    if (event_integrator%z(1) /= 1.0_dp) then
      error stop 'boundary event endpoint is not on the boundary'
    end if
    if (.not. boundary_event_converged( &
        event_integrator%last_event_fraction_width, &
        event_integrator%last_event_radial_residual, 1.0e-11_dp, &
        1.0e-10_dp)) then
      error stop 'reported boundary event does not meet its tolerances'
    end if
    event_time = (real(step - 1, dp) + &
      event_integrator%last_step_fraction)*timestep
  end subroutine locate_linear_boundary
end program test_newton_solver_status
