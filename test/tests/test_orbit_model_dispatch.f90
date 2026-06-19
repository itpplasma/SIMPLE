program test_orbit_model_dispatch
  ! Wave-1 followup #1: orbit_model is read from the config namelist and the
  ! macrostep dispatch maps it to the right pusher. This test writes a minimal
  ! namelist, parses it via params%read_config, and asserts:
  !   - orbit_model is parsed (default 0 = GC; here set to ORBIT_PAULI).
  !   - cpp_stages_from_mode maps GAUSS1..4 to stage counts 1..4 (the CPP branch
  !     dispatch key), proving the integer-coded select-case path is wired.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use params, only: orbit_model, integmode, read_config
  use orbit_full, only: ORBIT_GC, ORBIT_PAULI, ORBIT_BORIS, ORBIT_FOSYMPL
  use orbit_symplectic_base, only: GAUSS1, GAUSS2, GAUSS3, GAUSS4
  use orbit_cpp, only: cpp_stages_from_mode

  implicit none

  integer :: nfail, u
  character(256) :: cfgfile

  nfail = 0
  cfgfile = 'test_orbit_model_dispatch.in'

  open (newunit=u, file=cfgfile, status='replace', action='write')
  write (u, '(A)') '&config'
  write (u, '(A)') '  orbit_model = 1'
  write (u, '(A)') '  integmode = 4'
  write (u, '(A)') '/'
  close (u)

  call read_config(cfgfile)

  call check('orbit_model parsed as ORBIT_PAULI', orbit_model == ORBIT_PAULI, nfail)
  call check('integmode parsed as GAUSS1', integmode == GAUSS1, nfail)

  ! The dispatch keys are distinct integers (no overlap).
  call check('orbit model codes distinct', &
      ORBIT_GC == 0 .and. ORBIT_PAULI == 1 .and. ORBIT_BORIS == 2 .and. &
      ORBIT_FOSYMPL == 3, nfail)

  ! Stage mapping that the CPP select-case dispatch uses.
  call check('GAUSS1 -> 1 stage', cpp_stages_from_mode(GAUSS1) == 1, nfail)
  call check('GAUSS2 -> 2 stages', cpp_stages_from_mode(GAUSS2) == 2, nfail)
  call check('GAUSS3 -> 3 stages', cpp_stages_from_mode(GAUSS3) == 3, nfail)
  call check('GAUSS4 -> 4 stages', cpp_stages_from_mode(GAUSS4) == 4, nfail)

  if (nfail == 0) then
    print *, 'ALL ORBIT-MODEL DISPATCH TESTS PASSED'
  else
    print *, 'ORBIT-MODEL DISPATCH TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine check(name, ok, nfail)
    character(*), intent(in) :: name
    logical, intent(in) :: ok
    integer, intent(inout) :: nfail
    if (ok) then
      print '(A,A)', 'PASS  ', name
    else
      print '(A,A)', 'FAIL  ', name
      nfail = nfail + 1
    end if
  end subroutine check

end program test_orbit_model_dispatch
