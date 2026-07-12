program test_newton_solver_status
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use field_can_mod, only: field_can_t
  use orbit_symplectic, only: guard_lobatto_stage_radii, newton_midpoint, &
    solve_newton_system
  use orbit_symplectic_base, only: symplectic_integrator_t, &
    SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_OUTSIDE_DOMAIN, &
    SYMPLECTIC_STEP_MAXITER, SYMPLECTIC_STEP_LINEAR_SOLVE

  implicit none

  type(symplectic_integrator_t) :: integrator
  type(field_can_t) :: field
  real(dp) :: x(5), xlast(5), matrix(2, 2), residual(2)
  real(dp) :: lobatto_state(10), expected_lobatto_state(10)
  integer :: i, status

  field%Aph = 0.0_dp
  field%ro0 = 1.0_dp
  x = [0.5_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.5_dp]
  call newton_midpoint(integrator, field, x, 1.0e-15_dp, 1.0e-12_dp, &
    0, xlast, status)
  if (status /= SYMPLECTIC_STEP_MAXITER) then
    error stop 'zero-iteration midpoint solve did not report max iterations'
  end if
  if (any(xlast /= x)) error stop 'zero-iteration solve changed the iterate'

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
end program test_newton_solver_status
