program test_cpp_invariants
  ! REFACTOR / CODE-MOTION ORACLE -- NOT a physics cross-validation.
  !
  ! Invariant conservation for the flux-canonical CPP pusher (orbit_cpp) on the
  ! BOOZER chart of the QA wout. This CPP is the GC degenerate-Lagrangian scheme
  ! with mu held fixed, so its residual is byte-identical to GC and its
  ! conserved-quantity behavior matches the validated GC integrator by
  ! construction. The "<= GC" asserts below are therefore IDENTITIES that guard
  ! the device-portable Newton/LU realization, not evidence that two distinct
  ! methods conserve equally well. The genuine 6D-Pauli invariant test (real
  ! gyration, mu adiabatically conserved, energy bounded) is
  ! test_cpp_pauli_gc_banana.
  !
  ! Asserts:
  !   - mu is a fixed parameter -> identically conserved (byte ==).
  !   - energy H = vpar^2/2 + mu*Bmod oscillation, energy secular drift, and the
  !     canonical p_phi excursion each match GC to a tight relative margin
  !     (identity guard), and are bounded with no secular energy growth.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: twopi
  use simple_main, only: init_field
  use simple, only: tracer_t, init_sympl, init_params
  use params, only: isw_field_type, field_input, coord_input, integmode
  use new_vmec_stuff_mod, only: rmajor
  use magfie_sub, only: BOOZER
  use field_can_mod, only: field_can_t, evaluate, get_val
  use orbit_symplectic_base, only: symplectic_integrator_t, GAUSS1
  use orbit_symplectic, only: orbit_timestep_sympl
  use orbit_cpp, only: orbit_timestep_cpp, orbit_cpp_init, cpp_stages_from_mode

  implicit none

  integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
  integer :: nfail
  type(tracer_t) :: norb

  nfail = 0

  isw_field_type = BOOZER
  field_input = 'wout.nc'
  coord_input = 'wout.nc'
  call init_field(norb, 'wout.nc', ans_s, ans_tp, amultharm, GAUSS1)
  call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0e-12_dp)
  if (.not. associated(evaluate)) then
    print *, 'evaluate pointer not associated for BOOZER field'
    error stop 1
  end if

  call run_invariants(norb, nfail)

  if (nfail == 0) then
    print *, 'ALL CPP-FLUX CODE-MOTION INVARIANT ORACLE TESTS PASSED'
  else
    print *, 'CPP-FLUX CODE-MOTION INVARIANT ORACLE TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine run_invariants(norb, nfail)
    type(tracer_t), intent(inout) :: norb
    integer, intent(inout) :: nfail

    real(dp) :: z0(5), dtau, relerr, rbig
    real(dp) :: mu_cpp, mu_cpp_end
    real(dp) :: osc_cpp, sec_cpp, pdev_cpp
    real(dp) :: osc_gc,  sec_gc,  pdev_gc
    integer :: s, nstep

    rbig = rmajor*1.0e2_dp
    dtau    = twopi*rbig/256.0_dp
    relerr  = 1.0e-13_dp
    nstep   = 3000
    s = cpp_stages_from_mode(GAUSS1)

    z0 = [0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.1_dp]
    integmode = GAUSS1
    norb%relerr = relerr

    call run_one(z0, dtau, relerr, s, .true.,  osc_cpp, sec_cpp, pdev_cpp, &
                 mu_cpp, mu_cpp_end)
    call run_one(z0, dtau, relerr, s, .false., osc_gc,  sec_gc,  pdev_gc, &
                 mu_cpp, mu_cpp_end)

    print '(A,ES12.4)', '  mu (fixed) = ', mu_cpp
    print '(A,ES12.4,A,ES12.4)', '  energy osc: CPP=', osc_cpp, ' GC=', osc_gc
    print '(A,ES12.4,A,ES12.4)', '  energy secular: CPP=', sec_cpp, ' GC=', sec_gc
    print '(A,ES12.4,A,ES12.4)', '  pphi excursion: CPP=', pdev_cpp, ' GC=', pdev_gc

    ! mu is a fixed parameter: it must not move at all.
    call check('mu identically conserved', mu_cpp == mu_cpp_end, nfail)
    ! no secular energy growth (structure preservation).
    call check('energy no secular drift', sec_cpp < 0.5_dp*osc_cpp + 1.0e-12_dp, &
               nfail)
    ! CPP conserves no worse than the validated GC integrator (relative margin).
    call check('energy osc <= GC', osc_cpp <= osc_gc*(1.0_dp + 1.0e-6_dp) + 1.0e-12_dp, &
               nfail)
    call check('energy secular <= GC', sec_cpp <= sec_gc*(1.0_dp + 1.0e-6_dp) + 1.0e-12_dp, &
               nfail)
    call check('pphi excursion ~ GC', &
               abs(pdev_cpp - pdev_gc) <= 1.0e-6_dp*max(pdev_gc, 1.0e-12_dp) + 1.0e-10_dp, &
               nfail)
  end subroutine run_invariants

  subroutine run_one(z0, dtau, relerr, s, use_cpp, osc, sec, pdev, mu0, muend)
    real(dp), intent(in) :: z0(5), dtau, relerr
    integer, intent(in) :: s
    logical, intent(in) :: use_cpp
    real(dp), intent(out) :: osc, sec, pdev, mu0, muend

    type(symplectic_integrator_t) :: si
    type(field_can_t) :: f
    real(dp) :: H0, H, Hmin, Hmax, pphi0, H_first, H_last
    integer :: it, ierr, nhalf, n1, n2

    call init_sympl(si, f, z0, dtau, dtau, relerr, GAUSS1)
    if (use_cpp) call orbit_cpp_init(si, f)

    mu0 = f%mu
    pphi0 = si%z(4)
    call get_val(f, si%z(4))
    H0 = f%H; Hmin = H0; Hmax = H0; pdev = 0.0_dp
    nhalf = 3000/2
    H_first = 0.0_dp; H_last = 0.0_dp; n1 = 0; n2 = 0

    do it = 1, 3000
      if (use_cpp) then
        call orbit_timestep_cpp(si, f, s, ierr)
      else
        call orbit_timestep_sympl(si, f, ierr)
      end if
      if (ierr /= 0) then
        print '(A,L1,A,I0)', '  early exit use_cpp=', use_cpp, ' at step ', it
        exit
      end if
      call get_val(f, si%z(4))
      H = f%H
      Hmin = min(Hmin, H); Hmax = max(Hmax, H)
      pdev = max(pdev, abs(si%z(4) - pphi0))
      if (it <= nhalf) then
        H_first = H_first + H; n1 = n1 + 1
      else
        H_last = H_last + H; n2 = n2 + 1
      end if
    end do

    muend = f%mu
    osc = (Hmax - Hmin)/abs(H0)
    sec = abs(H_last/max(n2,1) - H_first/max(n1,1))/abs(H0)
  end subroutine run_one

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

end program test_cpp_invariants
