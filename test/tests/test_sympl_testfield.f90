program test_sympl_testfield
  use, intrinsic :: iso_fortran_env, only : dp => real64
  use simple_main, only : init_field
  use simple, only : tracer_t, init_sympl
  use params, only : isw_field_type, field_input, integmode
  use magfie_sub, only : TEST
  use field_can_mod, only : evaluate
  use orbit_symplectic, only : orbit_timestep_sympl

  implicit none

  type(tracer_t) :: norb
  character(*), parameter :: vmec_file = WOUT_FILE
  integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
  integer :: ierr

  ! Configure symplectic GC with TEST field to ensure initialization succeeds
  isw_field_type = TEST
  field_input = ''
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
  call orbit_timestep_sympl(norb%si, norb%f, ierr)
  if (ierr /= 0) then
    print *, 'orbit_timestep_sympl failed with ierr=', ierr
    stop 1
  end if

  print *, 'TEST field symplectic step succeeded'
end program test_sympl_testfield
