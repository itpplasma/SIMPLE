module failed_symplectic_step_backend
  use field_can_base, only: field_can_t
  use orbit_symplectic_base, only: symplectic_integrator_t

  implicit none

contains

  subroutine fail_symplectic_step(si, f, ierr)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    integer, intent(out) :: ierr

    f%Bmod = 8.0d0
    f%vpar = 3.0d0
    f%mu = 1.0d0
    ierr = 1
  end subroutine fail_symplectic_step

end module failed_symplectic_step_backend

program test_sympl_testfield
  use, intrinsic :: iso_fortran_env, only : dp => real64, int64
  use failed_symplectic_step_backend, only: fail_symplectic_step
  use simple_main, only : init_field, macrostep, macrostep_with_wall_check
  use simple, only : tracer_t, init_sympl
  use params, only : isw_field_type, field_input, coord_input, integmode, &
    swcoll, orbit_model, ORBIT_GC
  use magfie_sub, only : TEST
  use field_can_mod, only : evaluate
  use orbit_symplectic, only : orbit_timestep_sympl
  use orbit_symplectic_base, only: SYMPLECTIC_STEP_OK

  implicit none

  type(tracer_t) :: norb
  character(*), parameter :: vmec_file = WOUT_FILE
  integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
  integer :: ierr
  real(dp) :: initial_si_z(4)

  ! Configure symplectic GC with TEST field to ensure initialization succeeds
  isw_field_type = TEST
  field_input = vmec_file
  coord_input = vmec_file
  integmode = 1

  call init_field(norb, vmec_file, ans_s, ans_tp, amultharm, integmode)

  if (.not. associated(evaluate)) then
    print *, 'evaluate pointer not associated for TEST field'
    stop 1
  end if

  ! Smoke step: one symplectic timestep with test field
  norb%relerr = 1.0e-12_dp
  norb%dtaumin = 1.0e-3_dp
  norb%dtau = 1.0e-3_dp

  call init_sympl(norb%si, norb%f, (/0.2_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/), &
                  norb%dtau, norb%dtaumin, norb%relerr, 1)

  ierr = 0
  initial_si_z = norb%si%z
  call orbit_timestep_sympl(norb%si, norb%f, ierr)
  if (ierr == SYMPLECTIC_STEP_OK) then
    error stop 'unconverged Euler test-field step lost its failure status'
  end if
  if (any(norb%si%z /= initial_si_z)) then
    error stop 'unconverged Euler test-field step changed accepted position'
  end if

  call test_failed_step_preserves_state

  print *, 'TEST field symplectic step succeeded'

contains

  subroutine test_failed_step_preserves_state
    real(dp), parameter :: initial_state(5) = [0.2_dp, 1.0_dp, 2.0_dp, &
      1.0_dp, 0.25_dp]
    real(dp) :: z(5)
    real(dp) :: x_previous(3)
    integer(int64) :: kt
    integer :: step_error

    swcoll = .false.
    orbit_model = ORBIT_GC
    orbit_timestep_sympl => fail_symplectic_step
    z = initial_state
    kt = 0_int64

    call macrostep(norb, z, kt, step_error, 1)

    if (step_error /= 1) error stop 'failed step status was not preserved'
    if (kt /= 0_int64) error stop 'failed step advanced the time index'
    if (any(z /= initial_state)) then
      error stop 'failed symplectic step changed the accepted state'
    end if

    z = initial_state
    x_previous = 0.0_dp
    kt = 0_int64
    call macrostep_with_wall_check(norb, z, kt, step_error, 1, 1, x_previous)
    if (step_error /= 1) error stop 'wall path lost failed step status'
    if (kt /= 0_int64) error stop 'failed wall path advanced the time index'
    if (any(z /= initial_state)) then
      error stop 'failed wall path changed the accepted state'
    end if
  end subroutine test_failed_step_preserves_state

end program test_sympl_testfield
